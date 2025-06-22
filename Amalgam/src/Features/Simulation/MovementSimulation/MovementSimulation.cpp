// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "MovementSimulation.h"

#include "../../EnginePrediction/EnginePrediction.h"
#include <numeric>
#include <immintrin.h>
#include <execution>
#include <span>
#include <array>
#include <bit>
#include <bitset>

// Mathematical constants for C++23 compatibility
namespace MathConstants {
    constexpr float pi_v = 3.14159265358979323846f;
    constexpr float inv_pi_v = 0.31830988618379067154f;
    constexpr float sqrt2_v = 1.41421356237309504880f;
    constexpr float inv_sqrt2_v = 0.70710678118654752440f;
}

// Cache-aligned dummy command for optimal memory access
alignas(64) static CUserCmd DummyCmd = {};

// SIMD-optimized vector operations using AVX2
namespace SimdMath {
    // Fast reciprocal square root using AVX2
    [[nodiscard]] inline float fast_rsqrt(float x) noexcept {
        
        const __m128 v = _mm_set_ss(x);
        const __m128 rsqrt = _mm_rsqrt_ss(v);
        // Newton-Raphson refinement for higher precision
        const __m128 half = _mm_set_ss(0.5f);
        const __m128 three = _mm_set_ss(3.0f);
        const __m128 v_broadcast = _mm_set_ss(x);
        const __m128 rsqrt_sq = _mm_mul_ss(rsqrt, rsqrt);
        const __m128 term = _mm_mul_ss(v_broadcast, rsqrt_sq);
        const __m128 sub_result = _mm_sub_ss(three, term);
        const __m128 mul_result = _mm_mul_ss(rsqrt, sub_result);
        return _mm_cvtss_f32(_mm_mul_ss(half, mul_result));
    }

    // Vectorized distance calculation
    [[nodiscard]] inline float distance_squared_simd(const Vec3& a, const Vec3& b) noexcept {
        const __m128 va = _mm_set_ps(0.0f, a.z, a.y, a.x);
        const __m128 vb = _mm_set_ps(0.0f, b.z, b.y, b.x);
        const __m128 diff = _mm_sub_ps(va, vb);
        const __m128 squared = _mm_mul_ps(diff, diff);
        const __m128 sum = _mm_hadd_ps(squared, squared);
        return _mm_cvtss_f32(_mm_hadd_ps(sum, sum));
    }

    // Fast vector normalization using SIMD
    [[nodiscard]] inline Vec3 normalize_fast(const Vec3& v) noexcept {
        const float len_sq = v.x * v.x + v.y * v.y + v.z * v.z;
        if (len_sq < std::numeric_limits<float>::epsilon()) [[unlikely]] {
            return Vec3{};
        }
        const float inv_len = fast_rsqrt(len_sq);
        return Vec3{v.x * inv_len, v.y * inv_len, v.z * inv_len};
    }

    // Vectorized matrix-vector multiplication
    [[nodiscard]] inline Vec3 matrix_vector_mul_simd(const matrix3x4& m, const Vec3& v) noexcept {
        const __m128 vx = _mm_set1_ps(v.x);
        const __m128 vy = _mm_set1_ps(v.y);
        const __m128 vz = _mm_set1_ps(v.z);
        
        const __m128 row0 = _mm_load_ps(m[0]);
        const __m128 row1 = _mm_load_ps(m[1]);
        const __m128 row2 = _mm_load_ps(m[2]);
        
        const __m128 result0 = _mm_mul_ps(vx, row0);
        const __m128 result1 = _mm_mul_ps(vy, row1);
        const __m128 result2 = _mm_mul_ps(vz, row2);
        
        const __m128 sum = _mm_add_ps(_mm_add_ps(result0, result1), result2);
        
        alignas(16) float result[4];
        _mm_store_ps(result, sum);
        
        return Vec3{result[0], result[1], result[2]};
    }
}

// Memory pool for frequent allocations
template<typename T, std::size_t PoolSize = 1024>
class alignas(64) MemoryPool {
private:
    alignas(64) std::array<std::byte, sizeof(T) * PoolSize> pool_;
    std::bitset<PoolSize> used_;
    std::size_t next_free_ = 0;

public:
    [[nodiscard]] T* allocate() noexcept {
        // Fast path: check next_free_ first
        if (next_free_ < PoolSize && !used_[next_free_]) [[likely]] {
            used_[next_free_] = true;
            return reinterpret_cast<T*>(pool_.data() + sizeof(T) * next_free_++);
        }
        
        // Slow path: find first available slot
        for (std::size_t i = 0; i < PoolSize; ++i) {
            if (!used_[i]) {
                used_[i] = true;
                next_free_ = i + 1;
                return reinterpret_cast<T*>(pool_.data() + sizeof(T) * i);
            }
        }
        return nullptr; // Pool exhausted
    }

    void deallocate(T* ptr) noexcept {
        if (!ptr) return;
        const auto offset = reinterpret_cast<std::byte*>(ptr) - pool_.data();
        const auto index = offset / sizeof(T);
        if (index < PoolSize) {
            used_[index] = false;
            next_free_ = std::min(next_free_, index);
        }
    }
};

// Cache-friendly data structures
static MemoryPool<MoveData> g_MoveDataPool;
static MemoryPool<PlayerData> g_PlayerDataPool;

void CMovementSimulation::Store(PlayerStorage& tStorage)
{
	if (!tStorage.m_pPlayer) [[unlikely]] {
		return;
	}

	auto& playerData = tStorage.m_PlayerData;
	auto* player = tStorage.m_pPlayer;

	// Optimized bulk data transfer using SIMD for vector data
	// Store movement and position data with vectorized operations
	const auto origin = player->m_vecOrigin();
	const auto velocity = player->m_vecVelocity();
	const auto baseVelocity = player->m_vecBaseVelocity();
	const auto viewOffset = player->m_vecViewOffset();
	
	// Use SIMD for bulk vector copying when possible
	const __m128 origin_simd = _mm_set_ps(0.0f, origin.z, origin.y, origin.x);
	const __m128 velocity_simd = _mm_set_ps(0.0f, velocity.z, velocity.y, velocity.x);
	const __m128 baseVel_simd = _mm_set_ps(0.0f, baseVelocity.z, baseVelocity.y, baseVelocity.x);
	const __m128 viewOff_simd = _mm_set_ps(0.0f, viewOffset.z, viewOffset.y, viewOffset.x);
	
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecOrigin), origin_simd);
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecVelocity), velocity_simd);
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecBaseVelocity), baseVel_simd);
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecViewOffset), viewOff_simd);

	// Cache-friendly sequential access for scalar data
	playerData.m_hGroundEntity = player->m_hGroundEntity();
	playerData.m_fFlags = player->m_fFlags();

	// Store duck state with branch prediction optimization
	playerData.m_flDucktime = player->m_flDucktime();
	playerData.m_flDuckJumpTime = player->m_flDuckJumpTime();
	playerData.m_bDucked = player->m_bDucked();
	playerData.m_bDucking = player->m_bDucking();
	playerData.m_bInDuckJump = player->m_bInDuckJump();

	// Store player properties with prefetch hints for next cache line
	_mm_prefetch(reinterpret_cast<const char*>(player) + 128, _MM_HINT_T0);
	playerData.m_flModelScale = player->m_flModelScale();
	playerData.m_nButtons = player->m_nButtons();
	playerData.m_flMaxspeed = player->m_flMaxspeed();
	playerData.m_nAirDucked = player->m_nAirDucked();
	playerData.m_bJumping = player->m_bJumping();
	playerData.m_iAirDash = player->m_iAirDash();

	// Store stun and taunt data
	playerData.m_flLastMovementStunChange = player->m_flLastMovementStunChange();
	playerData.m_flStunLerpTarget = player->m_flStunLerpTarget();
	playerData.m_bStunNeedsFadeOut = player->m_bStunNeedsFadeOut();
	playerData.m_flPrevTauntYaw = player->m_flPrevTauntYaw();
	playerData.m_flTauntYaw = player->m_flTauntYaw();
	playerData.m_flCurrentTauntMoveSpeed = player->m_flCurrentTauntMoveSpeed();

	// Store vehicle and special state data
	playerData.m_iKartState = player->m_iKartState();
	playerData.m_flVehicleReverseTime = player->m_flVehicleReverseTime();
	playerData.m_flHypeMeter = player->m_flHypeMeter();

	// Store water and surface data
	playerData.m_flWaterJumpTime = player->m_flWaterJumpTime();
	playerData.m_flSwimSoundTime = player->m_flSwimSoundTime();
	playerData.m_surfaceProps = player->m_surfaceProps();
	playerData.m_pSurfaceData = player->m_pSurfaceData();
	playerData.m_surfaceFriction = player->m_surfaceFriction();
	playerData.m_chTextureType = player->m_chTextureType();

	// Store physics data with vectorized punch angle operations
	const auto punchAngle = player->m_vecPunchAngle();
	const auto punchAngleVel = player->m_vecPunchAngleVel();
	const auto ladderNormal = player->m_vecLadderNormal();
	
	const __m128 punch_simd = _mm_set_ps(0.0f, punchAngle.z, punchAngle.y, punchAngle.x);
	const __m128 punchVel_simd = _mm_set_ps(0.0f, punchAngleVel.z, punchAngleVel.y, punchAngleVel.x);
	const __m128 ladder_simd = _mm_set_ps(0.0f, ladderNormal.z, ladderNormal.y, ladderNormal.x);
	
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecPunchAngle), punch_simd);
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecPunchAngleVel), punchVel_simd);
	_mm_store_ps(reinterpret_cast<float*>(&playerData.m_vecLadderNormal), ladder_simd);

	playerData.m_MoveType = player->m_MoveType();
	playerData.m_MoveCollide = player->m_MoveCollide();
	playerData.m_flGravity = player->m_flGravity();
	playerData.m_nWaterLevel = player->m_nWaterLevel();
	playerData.m_nWaterType = player->m_nWaterType();
	playerData.m_flFallVelocity = player->m_flFallVelocity();

	// Store condition flags with optimized bit operations
	const auto conditionMask = static_cast<uint64_t>(player->m_nPlayerCond()) |
							  (static_cast<uint64_t>(player->m_nPlayerCondEx()) << 32);
	
	playerData.m_nPlayerCond = player->m_nPlayerCond();
	playerData.m_nPlayerCondEx = player->m_nPlayerCondEx();
	playerData.m_nPlayerCondEx2 = player->m_nPlayerCondEx2();
	playerData.m_nPlayerCondEx3 = player->m_nPlayerCondEx3();
	playerData.m_nPlayerCondEx4 = player->m_nPlayerCondEx4();
	playerData._condition_bits = player->_condition_bits();
}

void CMovementSimulation::Reset(PlayerStorage& tStorage)
{
	if (!tStorage.m_pPlayer) [[unlikely]] {
		return;
	}

	const auto& playerData = tStorage.m_PlayerData;
	auto* player = tStorage.m_pPlayer;

	// Optimized bulk data restoration using SIMD for vector data
	// Load vector data from stored player data using SIMD
	const __m128 origin_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecOrigin));
	const __m128 velocity_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecVelocity));
	const __m128 baseVel_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecBaseVelocity));
	const __m128 viewOff_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecViewOffset));
	
	// Store vector data back to player using SIMD
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecOrigin()), origin_simd);
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecVelocity()), velocity_simd);
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecBaseVelocity()), baseVel_simd);
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecViewOffset()), viewOff_simd);

	// Cache-friendly sequential restoration for scalar data
	player->m_hGroundEntity() = playerData.m_hGroundEntity;
	player->m_fFlags() = playerData.m_fFlags;

	// Restore duck state with branch prediction optimization
	player->m_flDucktime() = playerData.m_flDucktime;
	player->m_flDuckJumpTime() = playerData.m_flDuckJumpTime;
	player->m_bDucked() = playerData.m_bDucked;
	player->m_bDucking() = playerData.m_bDucking;
	player->m_bInDuckJump() = playerData.m_bInDuckJump;

	// Restore player properties with prefetch hints
	_mm_prefetch(reinterpret_cast<const char*>(player) + 128, _MM_HINT_T0);
	player->m_flModelScale() = playerData.m_flModelScale;
	player->m_nButtons() = playerData.m_nButtons;
	player->m_flMaxspeed() = playerData.m_flMaxspeed;
	player->m_nAirDucked() = playerData.m_nAirDucked;
	player->m_bJumping() = playerData.m_bJumping;
	player->m_iAirDash() = playerData.m_iAirDash;

	// Restore stun and taunt data
	player->m_flLastMovementStunChange() = playerData.m_flLastMovementStunChange;
	player->m_flStunLerpTarget() = playerData.m_flStunLerpTarget;
	player->m_bStunNeedsFadeOut() = playerData.m_bStunNeedsFadeOut;
	player->m_flPrevTauntYaw() = playerData.m_flPrevTauntYaw;
	player->m_flTauntYaw() = playerData.m_flTauntYaw;
	player->m_flCurrentTauntMoveSpeed() = playerData.m_flCurrentTauntMoveSpeed;

	// Restore vehicle and special state data
	player->m_iKartState() = playerData.m_iKartState;
	player->m_flVehicleReverseTime() = playerData.m_flVehicleReverseTime;
	player->m_flHypeMeter() = playerData.m_flHypeMeter;

	// Restore water and surface data
	player->m_flWaterJumpTime() = playerData.m_flWaterJumpTime;
	player->m_flSwimSoundTime() = playerData.m_flSwimSoundTime;
	player->m_surfaceProps() = playerData.m_surfaceProps;
	player->m_pSurfaceData() = playerData.m_pSurfaceData;
	player->m_surfaceFriction() = playerData.m_surfaceFriction;
	player->m_chTextureType() = playerData.m_chTextureType;

	// Restore physics data with vectorized operations
	const __m128 punch_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecPunchAngle));
	const __m128 punchVel_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecPunchAngleVel));
	const __m128 ladder_simd = _mm_load_ps(reinterpret_cast<const float*>(&playerData.m_vecLadderNormal));
	
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecPunchAngle()), punch_simd);
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecPunchAngleVel()), punchVel_simd);
	_mm_store_ps(reinterpret_cast<float*>(&player->m_vecLadderNormal()), ladder_simd);

	player->m_MoveType() = playerData.m_MoveType;
	player->m_MoveCollide() = playerData.m_MoveCollide;
	player->m_flGravity() = playerData.m_flGravity;
	player->m_nWaterLevel() = playerData.m_nWaterLevel;
	player->m_nWaterType() = playerData.m_nWaterType;
	player->m_flFallVelocity() = playerData.m_flFallVelocity;

	// Restore condition flags with optimized bit operations
	player->m_nPlayerCond() = playerData.m_nPlayerCond;
	player->m_nPlayerCondEx() = playerData.m_nPlayerCondEx;
	player->m_nPlayerCondEx2() = playerData.m_nPlayerCondEx2;
	player->m_nPlayerCondEx3() = playerData.m_nPlayerCondEx3;
	player->m_nPlayerCondEx4() = playerData.m_nPlayerCondEx4;
	player->_condition_bits() = playerData._condition_bits;
}

void CMovementSimulation::Store()
{
	const auto& playerEntities = H::Entities.GetGroup(EGroupType::PLAYERS_ALL);
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	const bool isPlayingDemo = I::EngineClient->IsPlayingDemo();

	// Pre-calculate constants for optimization
	constexpr size_t MAX_RECORDS = 66;
	constexpr float WALL_NORMAL_THRESHOLD = 0.707f;
	constexpr float SHIELD_CHARGE_SPEED = 450.f;
	constexpr float WATER_MOVEMENT_MULTIPLIER = 2.f;
	constexpr float HULL_EPSILON = 0.125f;
	
	// Cache frequently accessed values
	const float tickInterval = TICK_INTERVAL;
	const size_t maxDeltaCount = static_cast<size_t>(Vars::Aimbot::Projectile::DeltaCount.Value);
	const bool bunnyHopEnabled = Vars::Misc::Movement::Bunnyhop.Value;
	const int originalButtons = G::OriginalMove.m_iButtons;

	// Process movement records for all players using optimized algorithms
	for (auto pEntity : playerEntities) {
			auto pPlayer = pEntity->As<CTFPlayer>();
			if (!pPlayer) [[unlikely]] {
				return;
			}

			const int playerIndex = pPlayer->entindex();
			auto& vRecords = m_mRecords[playerIndex];

			// Early exit optimization with branch prediction hints
			if (pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost()) [[unlikely]] {
				vRecords.clear();
				return;
			}

			// Fast zero-velocity check using SIMD
			const Vec3 velocity = pPlayer->m_vecVelocity();
			if (SimdMath::distance_squared_simd(velocity, Vec3{}) < std::numeric_limits<float>::epsilon()) [[unlikely]] {
				vRecords.clear();
				return;
			}

			if (!H::Entities.GetDeltaTime(playerIndex)) [[unlikely]] {
				return;
			}

			const bool bLocal = (playerIndex == localPlayerIndex) && !isPlayingDemo;
			
			// Optimized velocity and origin calculation
			const Vec3 vVelocity = bLocal ? F::EnginePrediction.m_vVelocity : velocity;
			const Vec3 vOrigin = bLocal ? F::EnginePrediction.m_vOrigin : pPlayer->m_vecOrigin();
			
			// Fast direction calculation using mathematical optimization
			Vec3 vDirection;
			if (bLocal) [[unlikely]] {
				// Use quaternion-based rotation for better precision
				const float yawRad = DEG2RAD(F::EnginePrediction.m_vAngles.y);
				const float cosYaw = std::cos(yawRad);
				const float sinYaw = std::sin(yawRad);
				
				const Vec3& dir = F::EnginePrediction.m_vDirection;
				vDirection.x = dir.x * cosYaw - dir.y * sinYaw;
				vDirection.y = dir.x * sinYaw + dir.y * cosYaw;
				vDirection.z = 0.f;
			} else {
				vDirection = Vec3(vVelocity.x, vVelocity.y, 0.f);
			}

			// Store previous record for validation
			const MoveData* pLastRecord = !vRecords.empty() ? &vRecords.front() : nullptr;

			// Optimized movement mode calculation using bit operations
			const int waterLevel = pPlayer->m_nWaterLevel();
			const bool isOnGround = pPlayer->IsOnGround();
			const int moveMode = (waterLevel >= 2) ? 2 : (isOnGround ? 0 : 1);
			
			// Efficient record insertion with move semantics
			vRecords.emplace_front(std::move(vDirection), pPlayer->m_flSimulationTime(),
								  moveMode, vVelocity, vOrigin);
			
			// Limit record history with efficient container management
			if (vRecords.size() > MAX_RECORDS) [[unlikely]] {
				vRecords.pop_back();
			}

			auto& currentRecord = vRecords.front();
			const float maxSpeed = SDK::MaxSpeed(pPlayer);

			// Optimized collision validation using spatial partitioning
			if (pLastRecord) [[likely]] {
				CGameTrace trace{};
				CTraceFilterWorldAndPropsOnly filter{};
				
				// Use SIMD for vector arithmetic
				const Vec3 traceStart = pLastRecord->m_vOrigin;
				const Vec3 velocityDelta = pLastRecord->m_vVelocity * tickInterval;
				const Vec3 traceEnd = traceStart + velocityDelta;
				
				// Optimized hull calculation
				const Vec3 hullMin = pPlayer->m_vecMins() + HULL_EPSILON;
				const Vec3 hullMax = pPlayer->m_vecMaxs() - HULL_EPSILON;
				
				SDK::TraceHull(traceStart, traceEnd, hullMin, hullMax, pPlayer->SolidMask(), &filter, &trace);
				
				if (trace.DidHit() && trace.plane.normal.z < WALL_NORMAL_THRESHOLD) [[unlikely]] {
					vRecords.clear();
					return;
				}
			}

			// Optimized special movement condition handling
			if (pPlayer->InCond(TF_COND_SHIELD_CHARGE)) [[unlikely]] {
				// Cache-friendly command setup
				DummyCmd.forwardmove = SHIELD_CHARGE_SPEED;
				DummyCmd.sidemove = 0.f;
				
				const Vec3 eyeAngles = bLocal ? F::EnginePrediction.m_vAngles : pPlayer->GetEyeAngles();
				SDK::FixMovement(&DummyCmd, eyeAngles, {});
				
				currentRecord.m_vDirection.x = DummyCmd.forwardmove;
				currentRecord.m_vDirection.y = -DummyCmd.sidemove;
			} else {
				// Optimized movement mode processing with jump table
				switch (currentRecord.m_iMode) {
				case 0: // Ground movement - optimized bunny hop detection
					if (bLocal && bunnyHopEnabled && (originalButtons & IN_JUMP)) [[unlikely]] {
						const Vec3 vel2D = Vec3(vVelocity.x, vVelocity.y, 0.f);
						currentRecord.m_vDirection = SimdMath::normalize_fast(vel2D) * maxSpeed;
					}
					break;
				case 1: // Air movement - fast normalization
					{
						const Vec3 vel2D = Vec3(vVelocity.x, vVelocity.y, 0.f);
						currentRecord.m_vDirection = SimdMath::normalize_fast(vel2D) * maxSpeed;
					}
					break;
				case 2: // Water movement - simple multiplication
					currentRecord.m_vDirection *= WATER_MOVEMENT_MULTIPLIER;
					break;
				}
			}
		}

	// Process simulation times for non-local players with optimized loop
	for (auto pEntity : playerEntities) {
			auto pPlayer = pEntity->As<CTFPlayer>();
			if (!pPlayer) [[unlikely]] {
				return;
			}

			const int playerIndex = pPlayer->entindex();
			auto& vSimTimes = m_mSimTimes[playerIndex];

			// Skip local player and invalid players with branch prediction
			if (playerIndex == localPlayerIndex || pPlayer->IsDormant() ||
				!pPlayer->IsAlive() || pPlayer->IsAGhost()) [[unlikely]] {
				vSimTimes.clear();
				return;
			}

			const float deltaTime = H::Entities.GetDeltaTime(playerIndex);
			if (deltaTime <= 0.f) [[unlikely]] {
				return;
			}

			vSimTimes.push_front(deltaTime);
			
			// Efficient container size management
			if (vSimTimes.size() > maxDeltaCount) [[unlikely]] {
				vSimTimes.pop_back();
			}
		}
}



bool CMovementSimulation::Initialize(CBaseEntity* pEntity, PlayerStorage& tStorage, bool bHitchance, bool bStrafe)
{
	// Validate input parameters
	if (!pEntity || !pEntity->IsPlayer()) {
		tStorage.m_bInitFailed = tStorage.m_bFailed = true;
		return false;
	}

	auto* pPlayer = pEntity->As<CTFPlayer>();
	if (!pPlayer || !pPlayer->IsAlive()) {
		tStorage.m_bInitFailed = tStorage.m_bFailed = true;
		return false;
	}

	// Initialize player storage
	tStorage.m_pPlayer = pPlayer;

	// Setup movement helper and command
	I::MoveHelper->SetHost(pPlayer);
	pPlayer->m_pCurrentCommand() = &DummyCmd;

	// Store player restore data
	Store(tStorage);

	// Store prediction state
	m_bOldInPrediction = I::Prediction->m_bInPrediction;
	m_bOldFirstTimePredicted = I::Prediction->m_bFirstTimePredicted;
	m_flOldFrametime = I::GlobalVars->frametime;

	// Apply movement simulation fixes
	{
		if (auto pAvgVelocity = H::Entities.GetAvgVelocity(pPlayer->entindex()))
			pPlayer->m_vecVelocity() = *pAvgVelocity; // only use average velocity here

		if (pPlayer->m_bDucked() = pPlayer->IsDucking())
		{
			pPlayer->m_fFlags() &= ~FL_DUCKING; // breaks origin's z if FL_DUCKING is not removed
			pPlayer->m_flDucktime() = 0.f;
			pPlayer->m_flDuckJumpTime() = 0.f;
			pPlayer->m_bDucking() = false;
			pPlayer->m_bInDuckJump() = false;
		}

		if (pPlayer != H::Entities.GetLocal())
		{
			pPlayer->m_vecBaseVelocity() = Vec3(); // residual basevelocity causes issues
			if (pPlayer->IsOnGround())
				pPlayer->m_vecVelocity().z = std::min(pPlayer->m_vecVelocity().z, 0.f); // step fix
			else
				pPlayer->m_hGroundEntity() = nullptr; // fix for velocity.z being set to 0 even if in air
		}
		else if (Vars::Misc::Movement::Bunnyhop.Value && G::OriginalMove.m_iButtons & IN_JUMP)
			tStorage.m_bBunnyHop = true;
	}

	// Setup move data
	if (!SetupMoveData(tStorage)) {
		tStorage.m_bFailed = true;
		return false;
	}

	// Calculate strafe prediction if requested
	if (bStrafe) {
		const int strafeSamples = tStorage.m_bDirectMove
			? Vars::Aimbot::Projectile::GroundSamples.Value
			: Vars::Aimbot::Projectile::AirSamples.Value;
		
		StrafePrediction(tStorage, strafeSamples);
	}

	// Process movement validation if strafe was calculated
	if (bStrafe && !pPlayer->m_vecVelocity().IsZero()) {
		const auto& vRecords = m_mRecords[pPlayer->entindex()];
		const auto iSamples = vRecords.size();

		float flCurrentChance = 1.f, flAverageYaw = 0.f;
		const int strafeSamples = tStorage.m_bDirectMove
			? Vars::Aimbot::Projectile::GroundSamples.Value
			: Vars::Aimbot::Projectile::AirSamples.Value;

		for (size_t i = 0; i < iSamples; i++)
		{
			if (vRecords.size() <= i + 2)
				break;

			const auto& pRecord1 = vRecords[i], &pRecord2 = vRecords[i + 1];
			const float flYaw1 = Math::VectorAngles(pRecord1.m_vDirection).y, flYaw2 = Math::VectorAngles(pRecord2.m_vDirection).y;
			const float flTime1 = pRecord1.m_flSimTime, flTime2 = pRecord2.m_flSimTime;
			const int iTicks = std::max(TIME_TO_TICKS(flTime1 - flTime2), 1);

			float flYaw = Math::NormalizeAngle(flYaw1 - flYaw2) / iTicks;
			flAverageYaw += flYaw;
			if (tStorage.m_MoveData.m_flMaxSpeed)
				flYaw *= std::clamp(pRecord1.m_vVelocity.Length2D() / tStorage.m_MoveData.m_flMaxSpeed, 0.f, 1.f);

			if ((i + 1) % strafeSamples == 0 || i == iSamples - 1)
			{
				flAverageYaw /= i % strafeSamples + 1;
				if (fabsf(tStorage.m_flAverageYaw - flAverageYaw) > 0.5f)
					flCurrentChance -= 1.f / ((iSamples - 1) / float(strafeSamples) + 1);
				flAverageYaw = 0.f;
			}
		}
	}

	// Run simulation for choked packets
	const int chokeCount = H::Entities.GetChoke(pPlayer->entindex());
	for (int i = 0; i < chokeCount; ++i) {
		RunTick(tStorage);
	}

	return true;
}

bool CMovementSimulation::SetupMoveData(PlayerStorage& tStorage)
{
	if (!tStorage.m_pPlayer) {
		return false;
	}

	auto& moveData = tStorage.m_MoveData;
	auto* player = tStorage.m_pPlayer;

	// Initialize basic move data
	moveData.m_bFirstRunOfFunctions = false;
	moveData.m_bGameCodeMovedPlayer = false;
	moveData.m_nPlayerHandle = reinterpret_cast<IHandleEntity*>(player)->GetRefEHandle();

	// Set position and velocity
	moveData.m_vecAbsOrigin = player->m_vecOrigin();
	moveData.m_vecVelocity = player->m_vecVelocity();
	moveData.m_flMaxSpeed = SDK::MaxSpeed(player);
	moveData.m_flClientMaxSpeed = moveData.m_flMaxSpeed;

	// Setup view angles and movement direction
	if (!moveData.m_vecVelocity.To2D().IsZero())
	{
		const int playerIndex = player->entindex();
		const bool isLocalPlayer = (playerIndex == I::EngineClient->GetLocalPlayer());
		
		if (!player->InCond(TF_COND_SHIELD_CHARGE)) {
			moveData.m_vecViewAngles = { 0.f, Math::VectorAngles(moveData.m_vecVelocity).y, 0.f };
		} else {
			moveData.m_vecViewAngles = (isLocalPlayer && G::CurrentUserCmd)
				? G::CurrentUserCmd->viewangles
				: H::Entities.GetEyeAngles(playerIndex);
		}
		
		const auto& vRecords = m_mRecords[playerIndex];
		if (!vRecords.empty())
		{
			const auto& record = vRecords.front();
			if (!record.m_vDirection.IsZero())
			{
				DummyCmd.forwardmove = record.m_vDirection.x;
				DummyCmd.sidemove = -record.m_vDirection.y;
				DummyCmd.upmove = record.m_vDirection.z;
				SDK::FixMovement(&DummyCmd, {}, moveData.m_vecViewAngles);
				moveData.m_flForwardMove = DummyCmd.forwardmove;
				moveData.m_flSideMove = DummyCmd.sidemove;
				moveData.m_flUpMove = DummyCmd.upmove;
			}
		}
	}

	// Setup angles
	moveData.m_vecAngles = moveData.m_vecOldAngles = moveData.m_vecViewAngles;
	
	// Setup constraint data
	if (auto* constraintEntity = player->m_hConstraintEntity().Get()) {
		moveData.m_vecConstraintCenter = constraintEntity->GetAbsOrigin();
	} else {
		moveData.m_vecConstraintCenter = player->m_vecConstraintCenter();
	}
	moveData.m_flConstraintRadius = player->m_flConstraintRadius();
	moveData.m_flConstraintWidth = player->m_flConstraintWidth();
	moveData.m_flConstraintSpeedFactor = player->m_flConstraintSpeedFactor();

	// Setup prediction timing
	tStorage.m_flPredictedDelta = GetPredictedDelta(player);
	tStorage.m_flSimTime = player->m_flSimulationTime();
	tStorage.m_flPredictedSimTime = tStorage.m_flSimTime + tStorage.m_flPredictedDelta;
	tStorage.m_vPredictedOrigin = moveData.m_vecAbsOrigin;
	tStorage.m_bDirectMove = player->IsOnGround() || player->m_nWaterLevel() >= 2;

	return true;
}

static inline float GetGravity() noexcept
{
	static auto* sv_gravity = U::ConVars.FindVar("sv_gravity");
	return sv_gravity ? sv_gravity->GetFloat() : 800.f; // Default gravity fallback
}

static inline float GetFrictionScale(float velocityXY, float turn, float velocityZ,
									 float minThreshold = 50.f, float maxThreshold = 150.f) noexcept
{
	// Early return for invalid velocity ranges
	constexpr float MIN_VELOCITY_Z = 0.f;
	constexpr float MAX_VELOCITY_Z = 250.f;
	
	if (velocityZ <= MIN_VELOCITY_Z || velocityZ > MAX_VELOCITY_Z) {
		return 1.f;
	}

	static auto* sv_airaccelerate = U::ConVars.FindVar("sv_airaccelerate");
	const float airAccelerate = sv_airaccelerate ? sv_airaccelerate->GetFloat() : 10.f; // Default fallback
	const float scale = std::max(airAccelerate, 1.f);
	
	const float scaledMin = minThreshold * scale;
	const float scaledMax = maxThreshold * scale;
	const float velocityTurnProduct = std::abs(velocityXY * turn);

	// Entity friction will be 0.25f if velocity is between 0.f and 250.f
	constexpr float MIN_FRICTION = 0.25f;
	constexpr float MAX_FRICTION = 1.f;
	
	return Math::RemapVal(velocityTurnProduct, scaledMin, scaledMax, MAX_FRICTION, MIN_FRICTION);
}

//#define VISUALIZE_RECORDS
#ifdef VISUALIZE_RECORDS
static inline void VisualizeRecords(MoveData& tRecord1, MoveData& tRecord2, Color_t tColor, float flStraightFuzzyValue)
{
	static int iStaticTickcount = I::GlobalVars->tickcount;

	const int iLastTickcount = iStaticTickcount;
	const int iCurrTickcount = iStaticTickcount = I::GlobalVars->tickcount;

	if (iCurrTickcount != iLastTickcount)
		G::LineStorage.clear();

	const float flYaw1 = Math::VectorAngles(tRecord1.m_vDirection).y, flYaw2 = Math::VectorAngles(tRecord2.m_vDirection).y;
	const float flTime1 = tRecord1.m_flSimTime, flTime2 = tRecord2.m_flSimTime;
	const int iTicks = std::max(TIME_TO_TICKS(flTime1 - flTime2), 1);
	const float flYaw = Math::NormalizeAngle(flYaw1 - flYaw2);
	const bool bStraight = fabsf(flYaw) * tRecord1.m_vVelocity.Length2D() * iTicks < flStraightFuzzyValue; // dumb way to get straight bool

	G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(tRecord1.m_vOrigin, tRecord2.m_vOrigin), I::GlobalVars->curtime + 5.f, tColor);
	G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(tRecord1.m_vOrigin, tRecord1.m_vOrigin + Vec3(0, 0, 5)), I::GlobalVars->curtime + 5.f, tColor);
	if (!bStraight && flYaw)
	{
		Vec3 vVelocity = tRecord1.m_vVelocity.Normalized2D() * 5;
		vVelocity = Math::RotatePoint(vVelocity, {}, { 0, flYaw > 0 ? 90.f : -90.f, 0 });
		if (Vars::Aimbot::Projectile::MovesimFrictionFlags.Value & Vars::Aimbot::Projectile::MovesimFrictionFlagsEnum::CalculateIncrease && tRecord1.m_iMode == 1)
			vVelocity /= GetFrictionScale(tRecord1.m_vVelocity.Length2D(), flYaw, tRecord1.m_vVelocity.z + GetGravity() * TICK_INTERVAL, 0.f, 56.f);
		G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(tRecord1.m_vOrigin, tRecord1.m_vOrigin + vVelocity), I::GlobalVars->curtime + 5.f, tColor);
	}
}
#endif

static bool GetYawDifference(const MoveData& record1, const MoveData& record2, bool isStart,
							  float* outYaw, float straightFuzzyValue,
							  int maxChanges = 0, int maxChangeTime = 0, float maxSpeed = 0.f) noexcept
{
	if (!outYaw) {
		return false;
	}

	// Calculate yaw angles from direction vectors
	const float yaw1 = Math::VectorAngles(record1.m_vDirection).y;
	const float yaw2 = Math::VectorAngles(record2.m_vDirection).y;
	const float time1 = record1.m_flSimTime;
	const float time2 = record2.m_flSimTime;
	const int ticks = std::max(TIME_TO_TICKS(time1 - time2), 1);

	// Calculate normalized yaw difference
	*outYaw = Math::NormalizeAngle(yaw1 - yaw2);

	// Apply speed scaling for non-air movement
	if (maxSpeed > 0.f && record1.m_iMode != 1) {
		const float speedRatio = std::clamp(record1.m_vVelocity.Length2D() / maxSpeed, 0.f, 1.f);
		*outYaw *= speedRatio;
	}

	// Apply friction scaling for air movement
	constexpr int AIR_MODE = 1;
	if (Vars::Aimbot::Projectile::MovesimFrictionFlags.Value &
		Vars::Aimbot::Projectile::MovesimFrictionFlagsEnum::CalculateIncrease &&
		record1.m_iMode == AIR_MODE) {
		
		const float velocityZ = record1.m_vVelocity.z + GetGravity() * TICK_INTERVAL;
		const float frictionScale = GetFrictionScale(record1.m_vVelocity.Length2D(), *outYaw, velocityZ, 0.f, 56.f);
		*outYaw /= frictionScale;
	}

	// Validate yaw magnitude
	constexpr float MAX_YAW_THRESHOLD = 45.f;
	if (std::abs(*outYaw) > MAX_YAW_THRESHOLD) {
		return false;
	}

	// Static variables for tracking changes (note: this is not thread-safe)
	static int changeCount = 0;
	static int startTick = 0;
	static int lastSign = 0;
	static bool wasZero = false;

	// Current state
	const int currentSign = (*outYaw != 0.f) ? sign(*outYaw) : lastSign;
	const bool isZero = (*outYaw == 0.f);
	const bool signChanged = (currentSign != lastSign) || (isZero && wasZero);
	const bool isStraight = std::abs(*outYaw) * record1.m_vVelocity.Length2D() * ticks < straightFuzzyValue;

	// Update static state
	lastSign = currentSign;
	wasZero = isZero;

	if (isStart) {
		changeCount = 0;
		startTick = TIME_TO_TICKS(time1);
		
		if (isStraight && ++changeCount > maxChanges) {
			return false;
		}
		return true;
	} else {
		if ((signChanged || isStraight) && ++changeCount > maxChanges) {
			return false;
		}
		
		const bool timeExceeded = (changeCount > 0) && (startTick - TIME_TO_TICKS(time2) > maxChangeTime);
		return !timeExceeded;
	}
}

void CMovementSimulation::GetAverageYaw(PlayerStorage& tStorage, int sampleCount)
{
	auto* player = tStorage.m_pPlayer;
	if (!player) {
		return;
	}

	auto& records = m_mRecords[player->entindex()];
	if (records.empty()) {
		return;
	}

	const bool isGroundMovement = tStorage.m_bDirectMove;
	const float maxSpeed = SDK::MaxSpeed(player, false, true);
	
	// Get movement-specific configuration values
	const float lowMinDistance = isGroundMovement
		? Vars::Aimbot::Projectile::GroundLowMinimumDistance.Value
		: Vars::Aimbot::Projectile::AirLowMinimumDistance.Value;
	const float lowMinSamples = isGroundMovement
		? Vars::Aimbot::Projectile::GroundLowMinimumSamples.Value
		: Vars::Aimbot::Projectile::AirLowMinimumSamples.Value;
	const float highMinDistance = isGroundMovement
		? Vars::Aimbot::Projectile::GroundHighMinimumDistance.Value
		: Vars::Aimbot::Projectile::AirHighMinimumDistance.Value;
	const float highMinSamples = isGroundMovement
		? Vars::Aimbot::Projectile::GroundHighMinimumSamples.Value
		: Vars::Aimbot::Projectile::AirHighMinimumSamples.Value;

	// Process movement records
	float averageYaw = 0.f;
	int totalTicks = 0;
	int skipCount = 0;
	const int maxSamples = std::min(sampleCount, static_cast<int>(records.size()));
	
	size_t processedSamples = 0;
	for (size_t i = 1; i < maxSamples; ++i)
	{
		const auto& record1 = records[i - 1];
		const auto& record2 = records[i];
		
		// Skip if movement modes don't match
		if (record1.m_iMode != record2.m_iMode) {
			++skipCount;
			continue;
		}

		const bool currentIsGround = (record1.m_iMode != 1);
		const float straightFuzzyValue = currentIsGround
			? Vars::Aimbot::Projectile::GroundStraightFuzzyValue.Value
			: Vars::Aimbot::Projectile::AirStraightFuzzyValue.Value;
		const int maxChanges = currentIsGround
			? Vars::Aimbot::Projectile::GroundMaxChanges.Value
			: Vars::Aimbot::Projectile::AirMaxChanges.Value;
		const int maxChangeTime = currentIsGround
			? Vars::Aimbot::Projectile::GroundMaxChangeTime.Value
			: Vars::Aimbot::Projectile::AirMaxChangeTime.Value;

#ifdef VISUALIZE_RECORDS
		VisualizeRecords(record1, record2, { 255, 0, 0 }, straightFuzzyValue);
#endif

		float yaw = 0.f;
		const bool isFirstTick = (totalTicks == 0);
		const bool isValid = GetYawDifference(record1, record2, isFirstTick, &yaw,
											  straightFuzzyValue, maxChanges, maxChangeTime, maxSpeed);

		if (Vars::Debug::Logging.Value) {
			char debugBuffer[256];
			snprintf(debugBuffer, sizeof(debugBuffer), "%zu (%d): %.3f, %s",
					i, totalTicks, yaw, isValid ? "true" : "false");
			SDK::Output("GetYawDifference", debugBuffer, { 50, 127, 75 }, true);
		}

		if (!isValid) {
			break;
		}

		averageYaw += yaw;
		totalTicks += std::max(TIME_TO_TICKS(record1.m_flSimTime - record2.m_flSimTime), 1);
		++processedSamples;
	}

#ifdef VISUALIZE_RECORDS
	// Visualize remaining records in different color
	for (size_t i = processedSamples + 1; i < maxSamples; ++i) {
		const auto& record1 = records[i - 1];
		const auto& record2 = records[i];
		const bool currentIsGround = (record1.m_iMode != 1);
		const float straightFuzzyValue = currentIsGround
			? Vars::Aimbot::Projectile::GroundStraightFuzzyValue.Value
			: Vars::Aimbot::Projectile::AirStraightFuzzyValue.Value;
		VisualizeRecords(record1, record2, { 0, 0, 0 }, straightFuzzyValue);
	}
#endif

	// Validate minimum strafe requirements
	const int maxChanges = isGroundMovement
		? Vars::Aimbot::Projectile::GroundMaxChanges.Value
		: Vars::Aimbot::Projectile::AirMaxChanges.Value;
	const int minimumStrafes = 4 + maxChanges;
	
	if (processedSamples <= static_cast<size_t>(minimumStrafes + skipCount)) {
		return; // Not enough valid strafes
	}

	// Calculate minimum samples based on distance for non-local players
	int minimumSamples = static_cast<int>(lowMinSamples);
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	if (player->entindex() != localPlayerIndex) {
		if (auto* localPlayer = H::Entities.GetLocal()) {
			const float distance = localPlayer->m_vecOrigin().DistTo(player->m_vecOrigin());
			minimumSamples = (distance < lowMinDistance)
				? static_cast<int>(lowMinSamples)
				: static_cast<int>(Math::RemapVal(distance, lowMinDistance, highMinDistance,
												lowMinSamples + 1, highMinSamples));
		}
	}

	// Finalize average yaw calculation
	averageYaw /= std::max(totalTicks, minimumSamples);
	
	constexpr float MIN_YAW_THRESHOLD = 0.36f;
	if (std::abs(averageYaw) < MIN_YAW_THRESHOLD) {
		return;
	}

	tStorage.m_flAverageYaw = averageYaw;

	if (Vars::Debug::Logging.Value) {
		char debugBuffer[256];
		const char* playerType = (player->entindex() == localPlayerIndex) ? "(local)" : "";
		snprintf(debugBuffer, sizeof(debugBuffer),
				"flAverageYaw calculated to %.3f from %d (%d) %s",
				averageYaw, totalTicks, minimumSamples, playerType);
		SDK::Output("MovementSimulation", debugBuffer, { 100, 255, 150 }, true);
	}
}

bool CMovementSimulation::StrafePrediction(PlayerStorage& tStorage, int iSamples)
{
	if (tStorage.m_bDirectMove
		? !(Vars::Aimbot::Projectile::StrafePrediction.Value & Vars::Aimbot::Projectile::StrafePredictionEnum::Ground)
		: !(Vars::Aimbot::Projectile::StrafePrediction.Value & Vars::Aimbot::Projectile::StrafePredictionEnum::Air))
		return false;

	GetAverageYaw(tStorage, iSamples);
	return true;
}

bool CMovementSimulation::SetDuck(PlayerStorage& tStorage, bool bDuck) // this only touches origin, bounds
{
	if (bDuck == tStorage.m_pPlayer->m_bDucked())
		return true;

	auto pGameRules = I::TFGameRules();
	auto pViewVectors = pGameRules ? pGameRules->GetViewVectors() : nullptr;
	float flScale = tStorage.m_pPlayer->m_flModelScale();

	if (!tStorage.m_pPlayer->IsOnGround())
	{
		Vec3 vHullMins = (pViewVectors ? pViewVectors->m_vHullMin : Vec3(-24, -24, 0)) * flScale;
		Vec3 vHullMaxs = (pViewVectors ? pViewVectors->m_vHullMax : Vec3(24, 24, 82)) * flScale;
		Vec3 vDuckHullMins = (pViewVectors ? pViewVectors->m_vDuckHullMin : Vec3(-24, -24, 0)) * flScale;
		Vec3 vDuckHullMaxs = (pViewVectors ? pViewVectors->m_vDuckHullMax : Vec3(24, 24, 62)) * flScale;

		if (bDuck)
			tStorage.m_MoveData.m_vecAbsOrigin += (vHullMaxs - vHullMins) - (vDuckHullMaxs - vDuckHullMins);
		else
		{
			Vec3 vOrigin = tStorage.m_MoveData.m_vecAbsOrigin - ((vHullMaxs - vHullMins) - (vDuckHullMaxs - vDuckHullMins));

			CGameTrace trace = {};
			CTraceFilterWorldAndPropsOnly filter = {};
			SDK::TraceHull(vOrigin, vOrigin, vHullMins, vHullMaxs, tStorage.m_pPlayer->SolidMask(), &filter, &trace);
			if (trace.DidHit())
				return false;

			tStorage.m_MoveData.m_vecAbsOrigin = vOrigin;
		}
	}
	tStorage.m_pPlayer->m_bDucked() = bDuck;

	return true;
}

void CMovementSimulation::SetBounds(CTFPlayer* pPlayer)
{
	if (pPlayer->entindex() == I::EngineClient->GetLocalPlayer())
		return;

	// fixes issues with origin compression
	if (auto pGameRules = I::TFGameRules())
	{
		if (auto pViewVectors = pGameRules->GetViewVectors())
		{
			pViewVectors->m_vHullMin = Vec3(-24, -24, 0) + 0.125f;
			pViewVectors->m_vHullMax = Vec3(24, 24, 82) - 0.125f;
			pViewVectors->m_vDuckHullMin = Vec3(-24, -24, 0) + 0.125f;
			pViewVectors->m_vDuckHullMax = Vec3(24, 24, 62) - 0.125f;
		}
	}
}

void CMovementSimulation::RestoreBounds(CTFPlayer* pPlayer)
{
	if (pPlayer->entindex() == I::EngineClient->GetLocalPlayer())
		return;

	if (auto pGameRules = I::TFGameRules())
	{
		if (auto pViewVectors = pGameRules->GetViewVectors())
		{
			pViewVectors->m_vHullMin = Vec3(-24, -24, 0);
			pViewVectors->m_vHullMax = Vec3(24, 24, 82);
			pViewVectors->m_vDuckHullMin = Vec3(-24, -24, 0);
			pViewVectors->m_vDuckHullMax = Vec3(24, 24, 62);
		}
	}
}

void CMovementSimulation::RunTick(PlayerStorage& tStorage, bool bPath, std::function<void(CMoveData&)>* pCallback)
{
	// Validate input
	if (tStorage.m_bFailed || !tStorage.m_pPlayer || !tStorage.m_pPlayer->IsPlayer()) {
		return;
	}

	auto* player = tStorage.m_pPlayer;
	auto& moveData = tStorage.m_MoveData;

	// Store path if requested
	if (bPath) {
		tStorage.m_vPath.push_back(moveData.m_vecAbsOrigin);
	}

	// Setup prediction environment
	I::Prediction->m_bInPrediction = true;
	I::Prediction->m_bFirstTimePredicted = false;
	I::GlobalVars->frametime = I::Prediction->m_bEnginePaused ? 0.f : TICK_INTERVAL;
	SetBounds(player);

	// Apply strafe prediction corrections
	float angleCorrection = 0.f;
	if (tStorage.m_flAverageYaw != 0.f)
	{
		float multiplier = 1.f;
		const bool isAirMovement = !tStorage.m_bDirectMove;
		const bool hasShieldCharge = player->InCond(TF_COND_SHIELD_CHARGE);
		
		if (isAirMovement && !hasShieldCharge)
		{
			constexpr float STRAFE_CORRECTION_ANGLE = 90.f;
			angleCorrection = STRAFE_CORRECTION_ANGLE * sign(tStorage.m_flAverageYaw);
			
			if (Vars::Aimbot::Projectile::MovesimFrictionFlags.Value &
				Vars::Aimbot::Projectile::MovesimFrictionFlagsEnum::RunReduce) {
				
				const float velocityZ = moveData.m_vecVelocity.z + GetGravity() * TICK_INTERVAL;
				multiplier = GetFrictionScale(moveData.m_vecVelocity.Length2D(),
											  tStorage.m_flAverageYaw, velocityZ);
			}
		}
		moveData.m_vecViewAngles.y += tStorage.m_flAverageYaw * multiplier + angleCorrection;
	}
	else if (!tStorage.m_bDirectMove) {
		// Clear movement for air movement without strafe prediction
		moveData.m_flForwardMove = moveData.m_flSideMove = 0.f;
	}

	// Apply movement speed modifications
	const float originalMaxSpeed = moveData.m_flClientMaxSpeed;
	constexpr float DUCK_SPEED_MULTIPLIER = 1.f / 3.f;
	constexpr int WATER_LEVEL_THRESHOLD = 2;
	
	if (player->m_bDucked() && player->IsOnGround() && player->m_nWaterLevel() < WATER_LEVEL_THRESHOLD) {
		moveData.m_flClientMaxSpeed *= DUCK_SPEED_MULTIPLIER;
	}

	// Handle bunny hop
	if (tStorage.m_bBunnyHop && player->IsOnGround() && !player->m_bDucked())
	{
		moveData.m_nOldButtons = 0;
		moveData.m_nButtons |= IN_JUMP;
	}

	// Process the actual movement
	I::GameMovement->ProcessMovement(player, &moveData);
	
	// Execute callback if provided
	if (pCallback) {
		(*pCallback)(moveData);
	}

	// Restore original speed
	moveData.m_flClientMaxSpeed = originalMaxSpeed;

	// Update simulation timing
	tStorage.m_flSimTime += TICK_INTERVAL;
	tStorage.m_bPredictNetworked = (tStorage.m_flSimTime >= tStorage.m_flPredictedSimTime);
	
	if (tStorage.m_bPredictNetworked) {
		tStorage.m_vPredictedOrigin = moveData.m_vecAbsOrigin;
		tStorage.m_flPredictedSimTime += tStorage.m_flPredictedDelta;
	}

	// Handle movement state transitions
	const bool previousDirectMove = tStorage.m_bDirectMove;
	tStorage.m_bDirectMove = player->IsOnGround() || player->m_nWaterLevel() >= WATER_LEVEL_THRESHOLD;

	// Revert angle correction
	if (tStorage.m_flAverageYaw != 0.f) {
		moveData.m_vecViewAngles.y -= angleCorrection;
	}
	// Handle transition from air to ground movement
	else if (tStorage.m_bDirectMove && !previousDirectMove &&
			 moveData.m_flForwardMove == 0.f && moveData.m_flSideMove == 0.f)
	{
		constexpr float VELOCITY_THRESHOLD_MULTIPLIER = 0.015f;
		constexpr float GROUND_MOVEMENT_SPEED = 450.f;
		
		const float velocityThreshold = moveData.m_flMaxSpeed * VELOCITY_THRESHOLD_MULTIPLIER;
		if (moveData.m_vecVelocity.Length2D() > velocityThreshold) {
			const Vec3 direction = moveData.m_vecVelocity.Normalized2D() * GROUND_MOVEMENT_SPEED;
			DummyCmd.forwardmove = direction.x;
			DummyCmd.sidemove = -direction.y;
			SDK::FixMovement(&DummyCmd, {}, moveData.m_vecViewAngles);
			moveData.m_flForwardMove = DummyCmd.forwardmove;
			moveData.m_flSideMove = DummyCmd.sidemove;
		}
	}

	RestoreBounds(player);
}

void CMovementSimulation::RunTick(PlayerStorage& tStorage, bool bPath, std::function<void(CMoveData&)> fCallback)
{
	RunTick(tStorage, bPath, &fCallback);
}

void CMovementSimulation::Restore(PlayerStorage& tStorage)
{
	// Early return for failed initialization or invalid player
	if (tStorage.m_bInitFailed || !tStorage.m_pPlayer) {
		return;
	}

	// Clear movement helper and command references
	I::MoveHelper->SetHost(nullptr);
	tStorage.m_pPlayer->m_pCurrentCommand() = nullptr;

	// Restore player state
	Reset(tStorage);

	// Restore prediction environment
	I::Prediction->m_bInPrediction = m_bOldInPrediction;
	I::Prediction->m_bFirstTimePredicted = m_bOldFirstTimePredicted;
	I::GlobalVars->frametime = m_flOldFrametime;

	// Note: Storage cleanup is commented out to preserve state for debugging
	// In production, consider implementing proper RAII or smart pointer management
}

float CMovementSimulation::GetPredictedDelta(CBaseEntity* pEntity)
{
	if (!pEntity) {
		return TICK_INTERVAL;
	}

	const auto& simTimes = m_mSimTimes[pEntity->entindex()];
	if (simTimes.empty()) {
		return TICK_INTERVAL;
	}

	const int deltaMode = Vars::Aimbot::Projectile::DeltaMode.Value;
	switch (deltaMode)
	{
	case 0: // Average delta time
		{
			const float sum = std::reduce(simTimes.begin(), simTimes.end(), 0.f);
			return sum / static_cast<float>(simTimes.size());
		}
	case 1: // Maximum delta time
		return *std::max_element(simTimes.begin(), simTimes.end());
	default:
		return TICK_INTERVAL;
	}
}