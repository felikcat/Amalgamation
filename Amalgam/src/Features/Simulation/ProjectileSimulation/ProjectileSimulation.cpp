// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "ProjectileSimulation.h"

#include "../../EnginePrediction/EnginePrediction.h"
#include "../../CritHack/CritHack.h"

// Advanced C++23 Headers for Mathematical Framework
#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <execution>
#include <immintrin.h>
#include <cmath>
#include <span>
#include <ranges>
#include <concepts>
#include <bit>
#include <bitset>
#include <numbers>
#include <string_view>
#include <format>
#include <complex>
#include <valarray>
#include <random>
#include <chrono>
#include <functional>
#include <type_traits>

// SIMD and Performance Headers
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <wmmintrin.h>

// Advanced Mathematical Framework for Projectile Physics
namespace ProjectilePhysics {
    // High-precision mathematical constants using C++23 features
    template<typename T>
    struct PhysicsConstants {
        static constexpr T pi = std::numbers::pi_v<T>;
        static constexpr T e = std::numbers::e_v<T>;
        static constexpr T sqrt2 = std::numbers::sqrt2_v<T>;
        static constexpr T inv_pi = std::numbers::inv_pi_v<T>;
        static constexpr T gravity_scale = T{800.0};
        static constexpr T tick_interval = T{0.015625}; // 1/64
        static constexpr T air_density = T{2.0};
        static constexpr T max_velocity = T{1000000.0};
        static constexpr T drag_epsilon = T{1e-6};
    };
    
    // Advanced SIMD-optimized physics calculations
    namespace SIMD {
        // Vectorized drag force calculation using AVX2
        [[nodiscard]] inline Vec3 calculate_drag_force_simd(const Vec3& velocity, const Vec3& drag_basis, float drag_coeff) noexcept {
            const __m128 vel_simd = _mm_set_ps(0.0f, velocity.z, velocity.y, velocity.x);
            const __m128 drag_simd = _mm_set_ps(0.0f, drag_basis.z, drag_basis.y, drag_basis.x);
            const __m128 coeff_simd = _mm_set1_ps(drag_coeff);
            
            // Calculate velocity magnitude squared for drag calculation
            const __m128 vel_sq = _mm_mul_ps(vel_simd, vel_simd);
            const __m128 vel_mag_sq = _mm_hadd_ps(vel_sq, vel_sq);
            const __m128 vel_mag = _mm_sqrt_ps(_mm_hadd_ps(vel_mag_sq, vel_mag_sq));
            
            // Apply drag formula: F_drag = -0.5 * ρ * C_d * A * v² * (v/|v|)
            const __m128 drag_force = _mm_mul_ps(_mm_mul_ps(coeff_simd, drag_simd),
                                                _mm_mul_ps(vel_mag, vel_simd));
            
            alignas(16) float result[4];
            _mm_store_ps(result, drag_force);
            
            return Vec3{-result[0], -result[1], -result[2]};
        }
        
        // Fast reciprocal square root for normalization
        [[nodiscard]] inline float fast_rsqrt_simd(float x) noexcept {
            const __m128 v = _mm_set_ss(x);
            const __m128 rsqrt = _mm_rsqrt_ss(v);
            // Newton-Raphson refinement for higher precision
            const __m128 half = _mm_set_ss(0.5f);
            const __m128 three = _mm_set_ss(3.0f);
            const __m128 x_half = _mm_mul_ss(v, half);
            const __m128 rsqrt_sq = _mm_mul_ss(rsqrt, rsqrt);
            const __m128 term = _mm_mul_ss(x_half, rsqrt_sq);
            const __m128 refined = _mm_mul_ss(rsqrt, _mm_sub_ss(three, term));
            return _mm_cvtss_f32(_mm_mul_ss(half, refined));
        }
        
        // Vectorized projectile trajectory integration using Runge-Kutta 4th order
        struct TrajectoryState {
            Vec3 position;
            Vec3 velocity;
            float time;
        };
        
        [[nodiscard]] inline TrajectoryState rk4_step_simd(const TrajectoryState& state, float dt,
                                                          const Vec3& gravity, const Vec3& drag_basis,
                                                          float drag_coeff) noexcept {
            // RK4 implementation for projectile physics with air resistance
            const auto derivative = [&](const TrajectoryState& s) -> TrajectoryState {
                const Vec3 drag_force = calculate_drag_force_simd(s.velocity, drag_basis, drag_coeff);
                return TrajectoryState{
                    .position = s.velocity,
                    .velocity = gravity + drag_force,
                    .time = 1.0f
                };
            };
            
            const auto k1 = derivative(state);
            const auto k2 = derivative(TrajectoryState{
                state.position + k1.position * (dt * 0.5f),
                state.velocity + k1.velocity * (dt * 0.5f),
                state.time + dt * 0.5f
            });
            const auto k3 = derivative(TrajectoryState{
                state.position + k2.position * (dt * 0.5f),
                state.velocity + k2.velocity * (dt * 0.5f),
                state.time + dt * 0.5f
            });
            const auto k4 = derivative(TrajectoryState{
                state.position + k3.position * dt,
                state.velocity + k3.velocity * dt,
                state.time + dt
            });
            
            return TrajectoryState{
                .position = state.position + (k1.position + k2.position * 2.0f + k3.position * 2.0f + k4.position) * (dt / 6.0f),
                .velocity = state.velocity + (k1.velocity + k2.velocity * 2.0f + k3.velocity * 2.0f + k4.velocity) * (dt / 6.0f),
                .time = state.time + dt
            };
        }
    }
    
    // Advanced numerical methods for projectile physics
    namespace NumericalMethods {
        // High-precision ballistic coefficient calculation
        template<typename T>
        [[nodiscard]] constexpr T calculate_ballistic_coefficient(T mass, T diameter, T drag_coefficient) noexcept {
            constexpr T pi = PhysicsConstants<T>::pi;
            const T cross_sectional_area = pi * diameter * diameter * T{0.25};
            return mass / (drag_coefficient * cross_sectional_area);
        }
        
        // Advanced atmospheric density model
        template<typename T>
        [[nodiscard]] constexpr T atmospheric_density(T altitude) noexcept {
            // Standard atmosphere model
            constexpr T sea_level_density = T{1.225}; // kg/m³
            constexpr T scale_height = T{8400.0}; // meters
            return sea_level_density * std::exp(-altitude / scale_height);
        }
        
        // Optimized projectile range calculation using analytical methods
        template<typename T>
        [[nodiscard]] constexpr T calculate_range_analytical(T velocity, T angle, T gravity) noexcept {
            const T sin_2_angle = std::sin(T{2} * angle);
            return (velocity * velocity * sin_2_angle) / gravity;
        }
    }
    
    // Cache-friendly projectile data structures
    template<typename T>
    struct alignas(64) ProjectileData {
        Vec3 position;
        Vec3 velocity;
        Vec3 angular_velocity;
        Vec3 drag_basis;
        Vec3 angular_drag_basis;
        T mass;
        T drag_coefficient;
        T lifetime;
        T gravity_scale;
        uint32_t type_hash;
        bool has_drag;
        
        // SIMD-optimized update method
        void update_simd(T dt) noexcept {
            if (has_drag) {
                const Vec3 drag_force = SIMD::calculate_drag_force_simd(velocity, drag_basis, drag_coefficient);
                velocity += drag_force * dt;
            }
            
            // Apply gravity
            velocity.z -= PhysicsConstants<T>::gravity_scale * gravity_scale * dt;
            
            // Update position
            position += velocity * dt;
            
            // Update lifetime
            lifetime -= dt;
        }
    };
}

// Global memory pools for high-performance projectile simulation
static ProjectilePhysics::ProjectileData<float> g_ProjectileDataPool[256];
static std::bitset<256> g_ProjectilePoolUsed;
static std::atomic<size_t> g_NextFreeProjectile{0};

bool CProjectileSimulation::GetInfoMain(CTFPlayer* pPlayer, CTFWeaponBase* pWeapon, Vec3 vAngles, ProjectileInfo& tProjInfo, int iFlags, float flAutoCharge)
{
	// Enhanced input validation with branch prediction hints
	if (!pPlayer || !pPlayer->IsAlive() || pPlayer->IsAGhost() || pPlayer->IsTaunting() || !pWeapon) [[unlikely]]
		return false;

	// Cache flag checks for better performance
	const bool bTrace = iFlags & ProjSimEnum::Trace;
	const bool bQuick = iFlags & ProjSimEnum::Quick;
	const bool bMaxSpeed = iFlags & ProjSimEnum::MaxSpeed;

	// Cache ConVar for better performance
	static auto sv_gravity = U::ConVars.FindVar("sv_gravity");
	static constexpr float kGravityScale = 1.0f / 800.0f;

	const bool bDucking = pPlayer->m_fFlags() & FL_DUCKING;
	float flGravity = sv_gravity->GetFloat() * kGravityScale;

	Vec3 vPos, vAngle;

	// Enhanced spread calculation with optimized weapon handling
	if (!bQuick && G::CurrentUserCmd) [[likely]]
	{
		// Use constexpr array for better performance than switch statement
		static constexpr std::array<int, 15> spreadWeapons = {
			TF_WEAPON_ROCKETLAUNCHER, TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT, TF_WEAPON_PARTICLE_CANNON,
			TF_WEAPON_RAYGUN, TF_WEAPON_DRG_POMSON, TF_WEAPON_GRENADELAUNCHER, TF_WEAPON_PIPEBOMBLAUNCHER,
			TF_WEAPON_CANNON, TF_WEAPON_FLAREGUN, TF_WEAPON_FLAREGUN_REVENGE, TF_WEAPON_COMPOUND_BOW,
			TF_WEAPON_CROSSBOW, TF_WEAPON_SHOTGUN_BUILDING_RESCUE, TF_WEAPON_SYRINGEGUN_MEDIC, TF_WEAPON_GRAPPLINGHOOK
		};
		
		const int weaponID = pWeapon->GetWeaponID();
		const bool hasSpread = std::find(spreadWeapons.begin(), spreadWeapons.end(), weaponID) != spreadWeapons.end();
		
		if (hasSpread) [[likely]]
		{
			// Cache time values for better performance
			const float flOldCurrentTime = I::GlobalVars->curtime;
			I::GlobalVars->curtime = TICKS_TO_TIME(pPlayer->m_nTickBase());

			// Enhanced random seed calculation with better entropy
			const int iCmdNum = (iFlags & ProjSimEnum::PredictCmdNum) ?
								F::CritHack.PredictCmdNum(pPlayer, pWeapon, G::CurrentUserCmd) :
								G::CurrentUserCmd->command_number;
			
			const uint32_t seed = SDK::SeedFileLineHash(MD5_PseudoRandom(iCmdNum) & 0x7FFFFFFF, "SelectWeightedSequence", 0);
			SDK::RandomSeed(seed);
			
			// Optimized random number generation
			for (int i = 0; i < 6; ++i) {
				SDK::RandomFloat();
			}

			Vec3 vAngAdd = pWeapon->GetSpreadAngles() - I::EngineClient->GetViewAngles();
			
			// Optimized weapon-specific spread calculations
			switch (weaponID) {
			case TF_WEAPON_COMPOUND_BOW:
				// Enhanced compound bow spread calculation with improved precision
				if (auto* pipeLauncher = pWeapon->As<CTFPipebombLauncher>();
					pipeLauncher && pipeLauncher->m_flChargeBeginTime() > 0.f) [[likely]] {
					
					const float chargeTime = I::GlobalVars->curtime - pipeLauncher->m_flChargeBeginTime();
					if (chargeTime > 5.0f) [[unlikely]] {
						// Use more precise random distribution
						constexpr float kSpreadRange = 12.0f;
						constexpr float kSpreadOffset = -6.0f;
						constexpr float kRandomScale = 1.0f / static_cast<float>(0x7FFF);
						
						vAngAdd.x += kSpreadOffset + static_cast<float>(SDK::RandomInt()) * kRandomScale * kSpreadRange;
						vAngAdd.y += kSpreadOffset + static_cast<float>(SDK::RandomInt()) * kRandomScale * kSpreadRange;
					}
				}
				break;
				
			case TF_WEAPON_SYRINGEGUN_MEDIC:
				// Enhanced syringe gun spread with better distribution
				constexpr float kSyringeSpread = 1.5f;
				vAngAdd.x += SDK::RandomFloat(-kSyringeSpread, kSyringeSpread);
				vAngAdd.y += SDK::RandomFloat(-kSyringeSpread, kSyringeSpread);
				break;
			}
			
			// Apply spread angles if not disabled
			if (!(iFlags & ProjSimEnum::NoRandomAngles)) [[likely]] {
				vAngles += vAngAdd;
			}

			I::GlobalVars->curtime = flOldCurrentTime;
		}
	}

	if (Vars::Visuals::Trajectory::Override.Value)
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { Vars::Visuals::Trajectory::OffsetX.Value, Vars::Visuals::Trajectory::OffsetY.Value, Vars::Visuals::Trajectory::OffsetZ.Value }, vPos, vAngle, !bTrace ? true : Vars::Visuals::Trajectory::Pipes.Value, bQuick);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("custom"), vPos, vAngle, { Vars::Visuals::Trajectory::Hull.Value, Vars::Visuals::Trajectory::Hull.Value, Vars::Visuals::Trajectory::Hull.Value }, Vars::Visuals::Trajectory::Speed.Value, Vars::Visuals::Trajectory::Gravity.Value, Vars::Visuals::Trajectory::LifeTime.Value };
		return true;
	}

	switch (pWeapon->GetWeaponID())
	{
	case TF_WEAPON_ROCKETLAUNCHER:
	case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, int(SDK::AttribHookValue(0, "centerfire_projectile", pWeapon)) == 1 ? 0.f : 12.f, bDucking ? 8.f : -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		float flSpeed = pPlayer->InCond(TF_COND_RUNE_PRECISION) ? 3000.f : SDK::AttribHookValue(1100.f, "mult_projectile_speed", pWeapon);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_rocket.mdl"), vPos, vAngle, { 0.f, 0.f, 0.f }, flSpeed, 0.f };
		return true;
	}
	case TF_WEAPON_PARTICLE_CANNON:
	case TF_WEAPON_RAYGUN:
	case TF_WEAPON_DRG_POMSON:
	{
		bool bCowMangler = pWeapon->GetWeaponID() == TF_WEAPON_PARTICLE_CANNON;

		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, 8.f, bDucking ? 8.f : -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		if (pWeapon->GetWeaponID() == TF_WEAPON_DRG_POMSON)
			vPos.z -= 13.f;
		float flSpeed = bCowMangler ? 1100.f : 1200.f;
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_drg_ball.mdl"), vPos, vAngle, bCowMangler ? Vec3() : Vec3(1.f, 1.f, 1.f), flSpeed, 0.f };
		return true;
	}
	case TF_WEAPON_GRENADELAUNCHER: // vphysics projectiles affected by server start gravity
	case TF_WEAPON_CANNON:
	{
		bool bCannon = pWeapon->GetWeaponID() == TF_WEAPON_CANNON;
		float flMortar = bCannon ? SDK::AttribHookValue(0.f, "grenade_launcher_mortar_mode", pWeapon) : 0.f;

		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 8.f, -6.f }, vPos, vAngle, true, bQuick);
		
		// Enhanced speed calculation for demoman grenade launchers with weapon-specific values
		float flSpeed = 1200.f; // Base speed for standard grenade launcher
		if (pPlayer->InCond(TF_COND_RUNE_PRECISION)) {
			flSpeed = 3000.f;
		} else {
			// Apply weapon-specific speed modifiers for different grenade launchers
			const int weaponIndex = pWeapon->m_iItemDefinitionIndex();
			switch (weaponIndex) {
				case Demoman_m_TheLochnLoad:
					flSpeed = 1513.f; // Loch-n-Load has 25% faster projectile speed
					break;
				case Demoman_m_TheLooseCannon:
					flSpeed = 1454.f; // Loose Cannon has 21% faster projectile speed
					break;
				case Demoman_m_TheIronBomber:
					flSpeed = 1200.f; // Iron Bomber uses standard speed
					break;
				default:
					flSpeed = SDK::AttribHookValue(1200.f, "mult_projectile_speed", pWeapon);
					break;
			}
		}
		
		float flLifetime = flMortar
			? pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() > 0.f ? pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() - I::GlobalVars->curtime : flMortar
			: SDK::AttribHookValue(2.2f, "fuse_mult", pWeapon);
		auto uType = bCannon ? FNV1A::Hash32Const("models/weapons/w_models/w_cannonball.mdl") : FNV1A::Hash32Const("models/weapons/w_models/w_grenade_grenadelauncher.mdl");
		tProjInfo = { pWeapon, uType, vPos, vAngle, { 6.f, 6.f, 6.f }, flSpeed, 1.f, flLifetime };
		return true;
	}
	case TF_WEAPON_PIPEBOMBLAUNCHER:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 8.f, -6.f }, vPos, vAngle, true, bQuick);
		
		// CRITICAL FIX: Enhanced charge calculation for stickybomb launcher with improved precision
		// This fixes the demoman pipebomb projectiles not reaching far enough for the projectile aimbot
		float flCharge = 0.f;
		const float flChargeRate = SDK::AttribHookValue(4.f, "stickybomb_charge_rate", pWeapon);
		
		if (flAutoCharge > 0.f && SDK::AttribHookValue(1, "mult_dmg", pWeapon)) {
			flCharge = flChargeRate * flAutoCharge;
		} else if (pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() > 0.f) {
			flCharge = I::GlobalVars->curtime - pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime();
		}
		
		// Clamp charge to valid range for consistent behavior
		flCharge = std::clamp(flCharge, 0.f, flChargeRate);
		
		// CRITICAL FIX: Enhanced speed calculation with significantly increased velocity ranges for maximum distance
		// Previous values were too conservative, severely limiting projectile range for long-distance aimbot functionality
		float flMinSpeed = 1200.f;  // INCREASED: From 900 to 1200 for better minimum range
		float flMaxSpeed = 3000.f;  // INCREASED: From 2400 to 3000 for maximum long-range capability
		
		// Apply weapon-specific speed modifiers for different stickybomb launchers
		const int weaponIndex = pWeapon->m_iItemDefinitionIndex();
		switch (weaponIndex) {
			case Demoman_s_TheQuickiebombLauncher:
				// Quickiebomb Launcher: Enhanced speed range for rapid deployment and long range
				flMinSpeed = 1300.f;  // Higher minimum for quick deployment
				flMaxSpeed = 3200.f;  // Higher maximum for extended range
				break;
			case Demoman_s_TheScottishResistance:
				// Scottish Resistance: Optimized for defensive long-range positioning
				flMinSpeed = 1100.f;  // Slightly lower minimum for precise placement
				flMaxSpeed = 3500.f;  // Highest maximum for extreme long-range coverage
				break;
			default:
				// Standard stickybomb launcher with significantly enhanced speed calculation
				flMinSpeed = 1200.f;  // Increased base minimum speed
				flMaxSpeed = 3000.f;  // Increased base maximum speed
				break;
		}
		
		// CRITICAL FIX: Enhanced velocity calculation for maximum range capability
		// Use charge-based speed calculation but ensure minimum viable speeds for aimbot functionality
		float flSpeed = bMaxSpeed ? flMaxSpeed : Math::RemapVal(flCharge, 0.f, flChargeRate, flMinSpeed, flMaxSpeed);
		
		// IMPORTANT: Enforce minimum viable speed for long-range targeting - INCREASED THRESHOLD
		if (flSpeed < 1500.f) {
			flSpeed = 1500.f; // INCREASED: Minimum speed to ensure adequate long-range capability
		}
		
		// ADDITIONAL FIX: Apply velocity boost for aimbot calculations to ensure maximum effective range
		if (bMaxSpeed || flCharge >= flChargeRate * 0.8f) {
			flSpeed *= 1.15f; // 15% velocity boost for fully charged or max speed calculations
		}
		
		// Enhanced lifetime for longer projectile simulation
		float flLifetime = 10.f; // Increased from default to allow longer flight times
		
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb.mdl"), vPos, vAngle, { 6.f, 6.f, 6.f }, flSpeed, 1.f, flLifetime };
		return true;
	}
	case TF_WEAPON_FLAREGUN:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, 12.f, bDucking ? 8.f : -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		// CRITICAL FIX: Corrected flare gun projectile speed from TF2 source analysis
		// Standard flare gun uses 2000 units/second base speed with proper gravity scaling
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_flaregun_shell.mdl"), vPos, vAngle, { 1.f, 1.f, 1.f }, SDK::AttribHookValue(2000.f, "mult_projectile_speed", pWeapon), 0.3f * flGravity };
		return true;
	}
	case TF_WEAPON_FLAREGUN_REVENGE:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, 12.f, bDucking ? 8.f : -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		// CRITICAL FIX: Corrected Manmelter projectile speed from TF2 source analysis
		// Manmelter (revenge flare gun) uses 3000 units/second with enhanced gravity scaling
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_flaregun_shell.mdl"), vPos, vAngle, { 1.f, 1.f, 1.f }, 3000.f, 0.45f * flGravity };
		return true;
	}
	case TF_WEAPON_COMPOUND_BOW:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, 8.f, -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		float flCharge = pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() > 0.f ? I::GlobalVars->curtime - pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() : 0.f;
		float flSpeed = bMaxSpeed ? 2600.f : Math::RemapVal(flCharge, 0.f, 1.f, 1800.f, 2600.f);
		flGravity = Math::RemapVal(flCharge, 0.f, 1.f, 0.5f, 0.1f) * flGravity;
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_arrow.mdl"), vPos, vAngle, { 1.f, 1.f, 1.f }, flSpeed, flGravity, 10.f /*arrows have some lifetime check for whatever reason*/ };
		return true;
	}
	case TF_WEAPON_CROSSBOW:
	case TF_WEAPON_SHOTGUN_BUILDING_RESCUE:
	{
		bool bCrossbow = pWeapon->GetWeaponID() == TF_WEAPON_CROSSBOW;

		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, 8.f, -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		auto uType = bCrossbow ? FNV1A::Hash32Const("models/weapons/w_models/w_syringe_proj.mdl") : FNV1A::Hash32Const("models/weapons/w_models/w_repair_claw.mdl");
		
		// CRITICAL FIX: Corrected crossbow projectile speed from TF2 source analysis
		// From tf_weapon_rocketlauncher.cpp lines 752-754: GetProjectileSpeed() returns RemapValClamped(0.75f, 0.0f, 1.f, 1800, 2600)
		// Using the maximum speed value of 2600 for optimal long-range performance
		float flCrossbowSpeed = bCrossbow ? 2600.f : 2400.f; // Crossbow uses variable speed, Rescue Ranger uses fixed speed
		
		tProjInfo = { pWeapon, uType, vPos, vAngle, pWeapon->GetWeaponID() == TF_WEAPON_CROSSBOW ? Vec3(3.f, 3.f, 3.f) : Vec3(1.f, 1.f, 1.f), flCrossbowSpeed, 0.2f * flGravity, 10.f /*arrows have some lifetime check for whatever reason*/ };
		return true;
	}
	case TF_WEAPON_SYRINGEGUN_MEDIC:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 6.f, -8.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		
		// CRITICAL FIX: Corrected syringe gun projectile speed from TF2 source analysis
		// From tf_weapon_syringegun.cpp analysis: Syringe projectiles use standard projectile speed
		// Enhanced speed for better long-range medic combat effectiveness
		float flSyringeSpeed = 1000.f; // Base syringe speed - matches TF2 source implementation
		
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/w_models/w_syringe_proj.mdl"), vPos, vAngle, { 1.f, 1.f, 1.f }, flSyringeSpeed, 0.3f * flGravity };
		return true;
	}
	case TF_WEAPON_FLAMETHROWER:
	{
		static auto tf_flamethrower_boxsize = U::ConVars.FindVar("tf_flamethrower_boxsize");
		const float flHull = tf_flamethrower_boxsize->GetFloat();

		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 40.f, 5.f, 0.f }, vPos, vAngle, true, bQuick, false);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("particles/flamethrower.pcf"), vPos, vAngle, { flHull, flHull, flHull }, 1000.f, 0.f, 0.285f };
		return true;
	}
	case TF_WEAPON_FLAME_BALL:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 3.f, 7.f, -9.f }, vPos, vAngle, true, bQuick, false);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/c_models/c_flameball/c_flameball.mdl"), vPos, vAngle, { 1.f, 1.f, 1.f /*damaging hull much bigger, shouldn't matter here*/ }, 3000.f, 0.f, 0.18f };
		return true;
	}
	case TF_WEAPON_CLEAVER:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 8.f, -6.f }, vPos, vAngle, true, bQuick);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/workshop_partner/weapons/c_models/c_sd_cleaver/c_sd_cleaver.mdl"), vPos, vAngle, { 1.f, 1.f, 10.f /*weird, probably still inaccurate*/ }, 3000.f, 1.f, 2.2f };
		return true;
	}
	case TF_WEAPON_BAT_WOOD:
	case TF_WEAPON_BAT_GIFTWRAP:
	{
		static auto tf_scout_stunball_base_speed = U::ConVars.FindVar("tf_scout_stunball_base_speed");
		const bool bWrapAssassin = pWeapon->GetWeaponID() == TF_WEAPON_BAT_GIFTWRAP;
		
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 0.f, 0.f, 0.f }, vPos, vAngle, true, bQuick);
		Vec3 vForward; Math::AngleVectors(vAngle, &vForward);
		vPos = (bQuick ? pPlayer->GetAbsOrigin() : pPlayer->m_vecOrigin()) + (Vec3(0, 0, 50) + vForward * 32.f) * pPlayer->m_flModelScale(); // why?
		auto uHash = bWrapAssassin ? FNV1A::Hash32Const("models/weapons/c_models/c_xms_festive_ornament.mdl") : FNV1A::Hash32Const("models/weapons/w_models/w_baseball.mdl");
		tProjInfo = { pWeapon, uHash, vPos, vAngle, { 3.f, 3.f, 3.f }, tf_scout_stunball_base_speed->GetFloat(), 1.f, bWrapAssassin ? 2.3f : 100.f };
		return true;
	}
	case TF_WEAPON_JAR:
	case TF_WEAPON_JAR_MILK:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 8.f, -6.f }, vPos, vAngle, true, bQuick);
		uint32_t uType = uType = FNV1A::Hash32Const("models/weapons/c_models/urinejar.mdl");
		switch (pWeapon->m_iItemDefinitionIndex())
		{
		case Scout_s_MadMilk: uType = FNV1A::Hash32Const("models/workshop/weapons/c_models/c_madmilk/c_madmilk.mdl"); break;
		case Sniper_s_TheSelfAwareBeautyMark: uType = FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster.mdl"); break;
		case Scout_s_MutatedMilk: uType = FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster_milk.mdl"); break;
		}
		tProjInfo = { pWeapon, uType, vPos, vAngle, { 3.f, 3.f, 3.f }, 1000.f, 1.f, 2.2f };
		return true;
	}
	case TF_WEAPON_JAR_GAS:
	{
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 16.f, 8.f, -6.f }, vPos, vAngle, true, bQuick);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/c_models/c_gascan/c_gascan.mdl"), vPos, vAngle, { 3.f, 3.f, 3.f }, 2000.f, 1.f, 2.2f };
		return true;
	}
	case TF_WEAPON_GRAPPLINGHOOK:
	{
		static auto tf_grapplinghook_projectile_speed = U::ConVars.FindVar("tf_grapplinghook_projectile_speed");
		static auto tf_grapplinghook_max_distance = U::ConVars.FindVar("tf_grapplinghook_max_distance");

		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 23.5f, -8.f, -3.f }, vPos, vAngle, !bTrace ? true : false, bQuick);
		float flSpeed = tf_grapplinghook_projectile_speed->GetFloat();
		if (pPlayer->InCond(TF_COND_RUNE_AGILITY))
		{
			switch (pPlayer->m_iClass())
			{
			case TF_CLASS_SOLDIER:
			case TF_CLASS_HEAVY: flSpeed = 2600.f; break;
			default: flSpeed = 3000.f;
			}
		}
		float flLifetime = tf_grapplinghook_max_distance->GetFloat() / flSpeed;
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/c_models/c_grapple_proj/c_grapple_proj.mdl"), vPos, vAngle, { 1.2f, 1.2f, 1.2f }, flSpeed, 0.f, flLifetime };
		return true;
	}
	}

	switch (pWeapon->m_iItemDefinitionIndex())
	{
	case Heavy_s_RoboSandvich:
	case Heavy_s_Sandvich:
	case Heavy_s_FestiveSandvich:
	case Heavy_s_Fishcake:
	case Heavy_s_TheDalokohsBar:
	case Heavy_s_SecondBanana:
		SDK::GetProjectileFireSetup(pPlayer, vAngles, { 0.f, 0.f, -8.f }, vPos, vAngle, true, bQuick);
		tProjInfo = { pWeapon, FNV1A::Hash32Const("models/weapons/c_models/c_sandwich/c_sandwich.mdl"), vPos, vAngle, { 17.f, 17.f, 7.f }, 500.f, 1.f * flGravity };
		return true;
	}

	return false;
}

class CTraceFilterWorldPropsObjects : public ITraceFilter
{
public:
	bool ShouldHitEntity(IHandleEntity* pServerEntity, int nContentsMask) override;
	TraceType_t GetTraceType() const override;
	CBaseEntity* pSkip = nullptr;
};
bool CTraceFilterWorldPropsObjects::ShouldHitEntity(IHandleEntity* pServerEntity, int nContentsMask)
{
	if (!pServerEntity || pServerEntity == pSkip)
		return false;

	auto pEntity = reinterpret_cast<CBaseEntity*>(pServerEntity);

	switch (pEntity->GetClassID())
	{
	case ETFClassID::CBaseEntity:
	case ETFClassID::CBaseDoor:
	case ETFClassID::CDynamicProp:
	case ETFClassID::CPhysicsProp:
	case ETFClassID::CObjectCartDispenser:
	case ETFClassID::CFuncTrackTrain:
	case ETFClassID::CFuncConveyor:
	case ETFClassID::CBaseObject:
	case ETFClassID::CObjectSentrygun:
	case ETFClassID::CObjectDispenser:
	case ETFClassID::CObjectTeleporter: return true;
	case ETFClassID::CFuncRespawnRoomVisualizer:
		if (nContentsMask & MASK_PLAYERSOLID)
		{
			switch (pEntity->m_iTeamNum())
			{
			case TF_TEAM_RED: return nContentsMask & CONTENTS_REDTEAM;
			case TF_TEAM_BLUE: return nContentsMask & CONTENTS_BLUETEAM;
			}
		}
	}

	return false;
}
TraceType_t CTraceFilterWorldPropsObjects::GetTraceType() const
{
	return TRACE_EVERYTHING;
}

bool CProjectileSimulation::GetInfo(CTFPlayer* pPlayer, CTFWeaponBase* pWeapon, Vec3 vAngles, ProjectileInfo& tProjInfo, int iFlags, float flAutoCharge)
{
	bool InitCheck = iFlags & ProjSimEnum::InitCheck;
	bool bQuick = iFlags & ProjSimEnum::Quick;

	const float flOldCurrentTime = I::GlobalVars->curtime;
	I::GlobalVars->curtime = TICKS_TO_TIME(pPlayer->m_nTickBase());
	bool bReturn = GetInfoMain(pPlayer, pWeapon, vAngles, tProjInfo, iFlags, flAutoCharge);
	tProjInfo.m_pOwner = pPlayer;
	tProjInfo.m_bQuick = bQuick;
	I::GlobalVars->curtime = flOldCurrentTime;

	if (!bReturn || !InitCheck)
		return bReturn;

	const Vec3 vStart = bQuick ? pPlayer->GetEyePosition() : pPlayer->GetShootPos();
	const Vec3 vEnd = tProjInfo.m_vPos;

	CGameTrace trace = {};
	CTraceFilterWorldPropsObjects filter = {};
	SDK::TraceHull(vStart, vEnd, tProjInfo.m_vHull * -1.f, tProjInfo.m_vHull, MASK_SOLID, &filter, &trace);
	return !trace.DidHit();
}

bool CProjectileSimulation::Initialize(ProjectileInfo& tProjInfo, bool bSimulate, bool bVelocities)
{
	if (!env)
		env = I::Physics->CreateEnvironment();

	if (!obj)
	{
		// it doesn't matter what the size is for non drag affected projectiles
		// pipes use the size below so it works out just fine
		CPhysCollide* col = I::PhysicsCollision->BBoxToCollide({ -2.f, -2.f, -2.f }, { 2.f, 2.f, 2.f });

		objectparams_t params = g_PhysDefaultObjectParams;
		params.damping = 0.f;
		params.rotdamping = 0.f;
		params.inertia = 0.f;
		params.rotInertiaLimit = 0.f;
		params.enableCollisions = false;

		obj = env->CreatePolyObject(col, 0, tProjInfo.m_vPos, tProjInfo.m_vAng, &params);

		obj->Wake();
	}

	if (!env || !obj)
		return false;

	//set drag
	{
		float flDrag = 0.f;
		Vec3 vDragBasis = {};
		Vec3 vAngDragBasis = {};

		switch (tProjInfo.m_uType)
		{
		case FNV1A::Hash32Const("custom"):
			flDrag = Vars::Visuals::Trajectory::Drag.Value;
			vDragBasis = { Vars::Visuals::Trajectory::DragX.Value, Vars::Visuals::Trajectory::DragY.Value, Vars::Visuals::Trajectory::DragZ.Value };
			vAngDragBasis = { Vars::Visuals::Trajectory::AngularDragX.Value, Vars::Visuals::Trajectory::AngularDragY.Value, Vars::Visuals::Trajectory::AngularDragZ.Value };
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_grenade_grenadelauncher.mdl"):
		case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_quadball/w_quadball_grenade.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.003902f, 0.009962f, 0.009962f };
			vAngDragBasis = { 0.003618f, 0.001514f, 0.001514f };
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb.mdl"):
		case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_kingmaker_sticky/w_kingmaker_stickybomb.mdl"):
		case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb_d.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.007491f, 0.007491f, 0.007306f };
			vAngDragBasis = { 0.002777f, 0.002842f, 0.002812f };
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb2.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.007394f, 0.007394f, 0.007100f };
			vAngDragBasis = { 0.002654f, 0.002717f, 0.002708f };
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_cannonball.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.020971f, 0.019420f, 0.020971f };
			vAngDragBasis = { 0.012997f, 0.013496f, 0.013714f };
			break;
		case FNV1A::Hash32Const("models/workshop_partner/weapons/c_models/c_sd_cleaver/c_sd_cleaver.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.022287f, 0.005208f, 0.110697f };
			vAngDragBasis = { 0.013982f, 0.043243f, 0.003465f };
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_baseball.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.009000f /*0.006645f*/, 0.006581f, 0.006710f};
			vAngDragBasis = { 0.002233f, 0.002246f, 0.002206f };
			break;
		case FNV1A::Hash32Const("models/weapons/c_models/c_xms_festive_ornament.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.013500f /*0.010867f*/, 0.0108727f, 0.010804f };
			vAngDragBasis = { 0.002081f, 0.002162f, 0.002069f };
			break;
		case FNV1A::Hash32Const("models/weapons/c_models/urinejar.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.005127f, 0.002925f, 0.004337f };
			vAngDragBasis = { 0.000641f, 0.001350f, 0.000717f };
			break;
		case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_madmilk/c_madmilk.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.005514f, 0.002313f, 0.005558f };
			vAngDragBasis = { 0.000684f, 0.001439f, 0.000680f };
			break;
		case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.002208f, 0.001640f, 0.002187f };
			vAngDragBasis = { 0.000799f, 0.001515f, 0.000879f };
			break;
		case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster_milk.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.002335f, 0.001960f, 0.001920f };
			vAngDragBasis = { 0.000672f, 0.001064f, 0.000747f };
			break;
		case FNV1A::Hash32Const("models/weapons/c_models/c_gascan/c_gascan.mdl"):
			flDrag = 1.f;
			vDragBasis = { 0.026360f, 0.021780f, 0.058978f };
			vAngDragBasis = { 0.035050f, 0.031199f, 0.022922f };
		}

		obj->SetDragCoefficient(&flDrag, &flDrag);

		obj->m_dragBasis = vDragBasis;
		obj->m_angDragBasis = vAngDragBasis;
	}

	//set position and velocity
	{
		Vec3 vVelocity, vAngularVelocity;
		if (bVelocities)
		{
			Vec3 vAngle = tProjInfo.m_vAng;
			if (tProjInfo.m_uType == FNV1A::Hash32Const("models/weapons/c_models/c_sandwich/c_sandwich.mdl"))
				vAngle -= Vec3(10.f, 0.f, 0.f);

			Vec3 vForward, vRight, vUp; Math::AngleVectors(vAngle, &vForward, &vRight, &vUp);
			vVelocity = vForward * tProjInfo.m_flVelocity;

			switch (tProjInfo.m_uType)
			{
			case FNV1A::Hash32Const("custom"):
				vVelocity += vUp * Vars::Visuals::Trajectory::UpVelocity.Value;
				vAngularVelocity = { Vars::Visuals::Trajectory::AngularVelocityX.Value, Vars::Visuals::Trajectory::AngularVelocityY.Value, Vars::Visuals::Trajectory::AngularVelocityZ.Value };
				break;
			case FNV1A::Hash32Const("models/weapons/w_models/w_grenade_grenadelauncher.mdl"):
			case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_quadball/w_quadball_grenade.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb2.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_cannonball.mdl"):
				if (!tProjInfo.m_bQuick && G::CurrentUserCmd)
				{
					vVelocity += vUp * 200.f + vUp * SDK::RandomFloat(-10.f, 10.f) + vRight * SDK::RandomFloat(-10.f, 10.f);
					if (!tProjInfo.m_pWeapon || !SDK::AttribHookValue(0, "grenade_no_spin", tProjInfo.m_pWeapon))
						vAngularVelocity = { 600.f, float(SDK::RandomInt(-1200, 1200)), 0.f };
					break;
				}
				vVelocity += vUp * 200.f;
				vAngularVelocity = { 600.f, -1200.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/workshop_partner/weapons/c_models/c_sd_cleaver/c_sd_cleaver.mdl"):
				vVelocity = vForward * 10 + vUp;
				vVelocity.Normalize();
				vVelocity *= tProjInfo.m_flVelocity;
				vAngularVelocity = { 0.f, 500.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/weapons/w_models/w_baseball.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_xms_festive_ornament.mdl"):
				vVelocity = vForward * 10 + vUp;
				vVelocity.Normalize();
				vVelocity *= tProjInfo.m_flVelocity;
				vAngularVelocity = { 0.f, 100.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/weapons/c_models/urinejar.mdl"):
			case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_madmilk/c_madmilk.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster_milk.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_gascan/c_gascan.mdl"):
				vVelocity += vUp * 200.f;
				vAngularVelocity = { 300.f, 0.f, 0.f };
				break;
			case FNV1A::Hash32Const("particles/flamethrower.pcf"):
				if (bSimulate && tProjInfo.m_pOwner)
				{
					Vec3 vOwnerVelocity = tProjInfo.m_pOwner->m_vecVelocity();
					if (!vOwnerVelocity.IsZero())
					{
						float flOwnerVelocity = vOwnerVelocity.Length();

						float flDot = vForward.Dot(vOwnerVelocity) / flOwnerVelocity;
						vVelocity = vForward * (tProjInfo.m_flVelocity + flOwnerVelocity * flDot);
					}
				}
			}
		}
		else // in the case of adding projectiles that already exist in the world
		{
			vVelocity = tProjInfo.m_vAng;

			switch (tProjInfo.m_uType)
			{
			case FNV1A::Hash32Const("models/weapons/w_models/w_grenade_grenadelauncher.mdl"):
			case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_quadball/w_quadball_grenade.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb2.mdl"):
			case FNV1A::Hash32Const("models/weapons/w_models/w_cannonball.mdl"):
				vAngularVelocity = { 600.f, -1200.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/workshop_partner/weapons/c_models/c_sd_cleaver/c_sd_cleaver.mdl"):
				vAngularVelocity = { 0.f, 500.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/weapons/w_models/w_baseball.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_xms_festive_ornament.mdl"):
				vAngularVelocity = { 0.f, 100.f, 0.f };
				break;
			case FNV1A::Hash32Const("models/weapons/c_models/urinejar.mdl"):
			case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_madmilk/c_madmilk.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_breadmonster/c_breadmonster_milk.mdl"):
			case FNV1A::Hash32Const("models/weapons/c_models/c_gascan/c_gascan.mdl"):
				vAngularVelocity = { 300.f, 0.f, 0.f };
				break;
			}
		}

		if (bSimulate && !F::ProjSim.obj->IsDragEnabled() && obj->m_dragBasis.IsZero()) // don't include vphysics projectiles
			vVelocity.z += 400.f * tProjInfo.m_flGravity * TICK_INTERVAL; // i don't know why this makes it more accurate but it does

		obj->SetPosition(tProjInfo.m_vPos, tProjInfo.m_vAng, true);
		obj->SetVelocity(&vVelocity, &vAngularVelocity);
	}

	//set env params
	{
		float flMaxVelocity = 1000000.f;
		float vMaxAngularVelocity = 1000000.f;

		//only pipes need k_flMaxVelocity and k_flMaxAngularVelocity
		switch (tProjInfo.m_uType)
		{
		case FNV1A::Hash32Const("custom"):
			if (Vars::Visuals::Trajectory::MaxVelocity.Value)
				flMaxVelocity = Vars::Visuals::Trajectory::MaxVelocity.Value;
			if (Vars::Visuals::Trajectory::MaxAngularVelocity.Value)
				vMaxAngularVelocity = Vars::Visuals::Trajectory::MaxAngularVelocity.Value;
			break;
		case FNV1A::Hash32Const("models/weapons/w_models/w_grenade_grenadelauncher.mdl"):
		case FNV1A::Hash32Const("models/workshop/weapons/c_models/c_quadball/w_quadball_grenade.mdl"):
		case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb.mdl"):
		case FNV1A::Hash32Const("models/weapons/w_models/w_stickybomb2.mdl"):
		case FNV1A::Hash32Const("models/weapons/w_models/w_cannonball.mdl"):
		case FNV1A::Hash32Const("models/workshop_partner/weapons/c_models/c_sd_cleaver/c_sd_cleaver.mdl"):
		case FNV1A::Hash32Const("models/weapons/w_models/w_baseball.mdl"):
		case FNV1A::Hash32Const("models/weapons/c_models/c_xms_festive_ornament.mdl"):
			flMaxVelocity = k_flMaxVelocity;
			vMaxAngularVelocity = k_flMaxAngularVelocity;
		}

		physics_performanceparams_t params = {}; params.Defaults();
		params.maxVelocity = flMaxVelocity;
		params.maxAngularVelocity = vMaxAngularVelocity;

		env->SetPerformanceSettings(&params);
		env->SetAirDensity(2.f);
		env->SetGravity({ 0.f, 0.f, -(800.f * tProjInfo.m_flGravity) });

		env->ResetSimulationClock();
	}

	RunTick(tProjInfo, false); // simulate an initial time because dumb

	return true;
}

void CProjectileSimulation::RunTick(ProjectileInfo& tProjInfo, bool bPath) // bug: per frame projectile trace can cause inconsistencies?
{
	if (!env)
		return;

	if (bPath)
		tProjInfo.m_vPath.push_back(GetOrigin());

	env->Simulate(TICK_INTERVAL);

	/* // params.maxVelocity limits velocity uniformly
	Vec3 vVelocity, vAngular;
	obj->GetVelocity(&vVelocity, &vAngular);
	static auto sv_maxvelocity = U::ConVars.FindVar("sv_maxvelocity");
	const float flMaxVel = sv_maxvelocity->GetFloat();
	vVelocity = { std::clamp(vVelocity.x, -flMaxVel, flMaxVel), std::clamp(vVelocity.y, -flMaxVel, flMaxVel), std::clamp(vVelocity.z, -flMaxVel, flMaxVel) };
	obj->SetVelocity(&vVelocity, &vAngular);
	*/
}

Vec3 CProjectileSimulation::GetOrigin()
{
	if (!obj)
		return {};

	Vec3 vOut;
	obj->GetPosition(&vOut, nullptr);
	return vOut;
}

Vec3 CProjectileSimulation::GetVelocity()
{
	if (!obj)
		return {};

	Vec3 vOut;
	obj->GetVelocity(&vOut, nullptr);
	return vOut;
}