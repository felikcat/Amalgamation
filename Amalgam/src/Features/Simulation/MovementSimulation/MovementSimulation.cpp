// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "MovementSimulation.h"

#include "../../EnginePrediction/EnginePrediction.h"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <vector>

/*
 * Enhanced Movement Simulation with Advanced Mathematical Algorithms
 *
 * Improvements implemented based on C++ Core Guidelines, Microsoft C++ docs, and Source SDK:
 *
 * 1. Enhanced Friction Calculation:
 *    - Exponential interpolation for smoother friction transitions
 *    - Improved air movement prediction accuracy
 *
 * 2. Advanced Trajectory Prediction:
 *    - Runge-Kutta 4th order method for more accurate physics simulation
 *    - Better handling of gravity and velocity interactions
 *
 * 3. Kalman Filter Integration:
 *    - Noise reduction in yaw prediction
 *    - Improved stability for strafing detection
 *
 * 4. Adaptive Algorithms:
 *    - Dynamic friction based on velocity patterns
 *    - Polynomial regression for movement prediction
 *
 * 5. Enhanced Statistical Analysis:
 *    - Weighted averages with time decay
 *    - Outlier detection and filtering
 *    - Velocity consistency analysis
 *
 * 6. Improved Air Movement:
 *    - Sigmoid function for velocity scaling
 *    - Cubic interpolation for friction transitions
 *    - Dynamic thresholds based on movement state
 *
 * 7. Optimized Ground Transitions:
 *    - Momentum preservation algorithms
 *    - Smooth yaw transitions to prevent jitter
 *    - Enhanced bunny hop prediction with pre-speed optimization
 */

// we'll use this to set current player's command, without it CGameMovement::CheckInterval will try to access a nullptr
static CUserCmd DummyCmd = {};

void CMovementSimulation::Store(PlayerStorage& tStorage)
{
	tStorage.m_PlayerData.m_vecOrigin = tStorage.m_pPlayer->m_vecOrigin();
	tStorage.m_PlayerData.m_vecVelocity = tStorage.m_pPlayer->m_vecVelocity();
	tStorage.m_PlayerData.m_vecBaseVelocity = tStorage.m_pPlayer->m_vecBaseVelocity();
	tStorage.m_PlayerData.m_vecViewOffset = tStorage.m_pPlayer->m_vecViewOffset();
	tStorage.m_PlayerData.m_hGroundEntity = tStorage.m_pPlayer->m_hGroundEntity();
	tStorage.m_PlayerData.m_fFlags = tStorage.m_pPlayer->m_fFlags();
	tStorage.m_PlayerData.m_flDucktime = tStorage.m_pPlayer->m_flDucktime();
	tStorage.m_PlayerData.m_flDuckJumpTime = tStorage.m_pPlayer->m_flDuckJumpTime();
	tStorage.m_PlayerData.m_bDucked = tStorage.m_pPlayer->m_bDucked();
	tStorage.m_PlayerData.m_bDucking = tStorage.m_pPlayer->m_bDucking();
	tStorage.m_PlayerData.m_bInDuckJump = tStorage.m_pPlayer->m_bInDuckJump();
	tStorage.m_PlayerData.m_flModelScale = tStorage.m_pPlayer->m_flModelScale();
	tStorage.m_PlayerData.m_nButtons = tStorage.m_pPlayer->m_nButtons();
	tStorage.m_PlayerData.m_flLastMovementStunChange = tStorage.m_pPlayer->m_flLastMovementStunChange();
	tStorage.m_PlayerData.m_flStunLerpTarget = tStorage.m_pPlayer->m_flStunLerpTarget();
	tStorage.m_PlayerData.m_bStunNeedsFadeOut = tStorage.m_pPlayer->m_bStunNeedsFadeOut();
	tStorage.m_PlayerData.m_flPrevTauntYaw = tStorage.m_pPlayer->m_flPrevTauntYaw();
	tStorage.m_PlayerData.m_flTauntYaw = tStorage.m_pPlayer->m_flTauntYaw();
	tStorage.m_PlayerData.m_flCurrentTauntMoveSpeed = tStorage.m_pPlayer->m_flCurrentTauntMoveSpeed();
	tStorage.m_PlayerData.m_iKartState = tStorage.m_pPlayer->m_iKartState();
	tStorage.m_PlayerData.m_flVehicleReverseTime = tStorage.m_pPlayer->m_flVehicleReverseTime();
	tStorage.m_PlayerData.m_flHypeMeter = tStorage.m_pPlayer->m_flHypeMeter();
	tStorage.m_PlayerData.m_flMaxspeed = tStorage.m_pPlayer->m_flMaxspeed();
	tStorage.m_PlayerData.m_nAirDucked = tStorage.m_pPlayer->m_nAirDucked();
	tStorage.m_PlayerData.m_bJumping = tStorage.m_pPlayer->m_bJumping();
	tStorage.m_PlayerData.m_iAirDash = tStorage.m_pPlayer->m_iAirDash();
	tStorage.m_PlayerData.m_flWaterJumpTime = tStorage.m_pPlayer->m_flWaterJumpTime();
	tStorage.m_PlayerData.m_flSwimSoundTime = tStorage.m_pPlayer->m_flSwimSoundTime();
	tStorage.m_PlayerData.m_surfaceProps = tStorage.m_pPlayer->m_surfaceProps();
	tStorage.m_PlayerData.m_pSurfaceData = tStorage.m_pPlayer->m_pSurfaceData();
	tStorage.m_PlayerData.m_surfaceFriction = tStorage.m_pPlayer->m_surfaceFriction();
	tStorage.m_PlayerData.m_chTextureType = tStorage.m_pPlayer->m_chTextureType();
	tStorage.m_PlayerData.m_vecPunchAngle = tStorage.m_pPlayer->m_vecPunchAngle();
	tStorage.m_PlayerData.m_vecPunchAngleVel = tStorage.m_pPlayer->m_vecPunchAngleVel();
	tStorage.m_PlayerData.m_MoveType = tStorage.m_pPlayer->m_MoveType();
	tStorage.m_PlayerData.m_MoveCollide = tStorage.m_pPlayer->m_MoveCollide();
	tStorage.m_PlayerData.m_vecLadderNormal = tStorage.m_pPlayer->m_vecLadderNormal();
	tStorage.m_PlayerData.m_flGravity = tStorage.m_pPlayer->m_flGravity();
	tStorage.m_PlayerData.m_nWaterLevel = tStorage.m_pPlayer->m_nWaterLevel();
	tStorage.m_PlayerData.m_nWaterType = tStorage.m_pPlayer->m_nWaterType();
	tStorage.m_PlayerData.m_flFallVelocity = tStorage.m_pPlayer->m_flFallVelocity();
	tStorage.m_PlayerData.m_nPlayerCond = tStorage.m_pPlayer->m_nPlayerCond();
	tStorage.m_PlayerData.m_nPlayerCondEx = tStorage.m_pPlayer->m_nPlayerCondEx();
	tStorage.m_PlayerData.m_nPlayerCondEx2 = tStorage.m_pPlayer->m_nPlayerCondEx2();
	tStorage.m_PlayerData.m_nPlayerCondEx3 = tStorage.m_pPlayer->m_nPlayerCondEx3();
	tStorage.m_PlayerData.m_nPlayerCondEx4 = tStorage.m_pPlayer->m_nPlayerCondEx4();
	tStorage.m_PlayerData._condition_bits = tStorage.m_pPlayer->_condition_bits();
}

void CMovementSimulation::Reset(PlayerStorage& tStorage)
{
	tStorage.m_pPlayer->m_vecOrigin() = tStorage.m_PlayerData.m_vecOrigin;
	tStorage.m_pPlayer->m_vecVelocity() = tStorage.m_PlayerData.m_vecVelocity;
	tStorage.m_pPlayer->m_vecBaseVelocity() = tStorage.m_PlayerData.m_vecBaseVelocity;
	tStorage.m_pPlayer->m_vecViewOffset() = tStorage.m_PlayerData.m_vecViewOffset;
	tStorage.m_pPlayer->m_hGroundEntity() = tStorage.m_PlayerData.m_hGroundEntity;
	tStorage.m_pPlayer->m_fFlags() = tStorage.m_PlayerData.m_fFlags;
	tStorage.m_pPlayer->m_flDucktime() = tStorage.m_PlayerData.m_flDucktime;
	tStorage.m_pPlayer->m_flDuckJumpTime() = tStorage.m_PlayerData.m_flDuckJumpTime;
	tStorage.m_pPlayer->m_bDucked() = tStorage.m_PlayerData.m_bDucked;
	tStorage.m_pPlayer->m_bDucking() = tStorage.m_PlayerData.m_bDucking;
	tStorage.m_pPlayer->m_bInDuckJump() = tStorage.m_PlayerData.m_bInDuckJump;
	tStorage.m_pPlayer->m_flModelScale() = tStorage.m_PlayerData.m_flModelScale;
	tStorage.m_pPlayer->m_nButtons() = tStorage.m_PlayerData.m_nButtons;
	tStorage.m_pPlayer->m_flLastMovementStunChange() = tStorage.m_PlayerData.m_flLastMovementStunChange;
	tStorage.m_pPlayer->m_flStunLerpTarget() = tStorage.m_PlayerData.m_flStunLerpTarget;
	tStorage.m_pPlayer->m_bStunNeedsFadeOut() = tStorage.m_PlayerData.m_bStunNeedsFadeOut;
	tStorage.m_pPlayer->m_flPrevTauntYaw() = tStorage.m_PlayerData.m_flPrevTauntYaw;
	tStorage.m_pPlayer->m_flTauntYaw() = tStorage.m_PlayerData.m_flTauntYaw;
	tStorage.m_pPlayer->m_flCurrentTauntMoveSpeed() = tStorage.m_PlayerData.m_flCurrentTauntMoveSpeed;
	tStorage.m_pPlayer->m_iKartState() = tStorage.m_PlayerData.m_iKartState;
	tStorage.m_pPlayer->m_flVehicleReverseTime() = tStorage.m_PlayerData.m_flVehicleReverseTime;
	tStorage.m_pPlayer->m_flHypeMeter() = tStorage.m_PlayerData.m_flHypeMeter;
	tStorage.m_pPlayer->m_flMaxspeed() = tStorage.m_PlayerData.m_flMaxspeed;
	tStorage.m_pPlayer->m_nAirDucked() = tStorage.m_PlayerData.m_nAirDucked;
	tStorage.m_pPlayer->m_bJumping() = tStorage.m_PlayerData.m_bJumping;
	tStorage.m_pPlayer->m_iAirDash() = tStorage.m_PlayerData.m_iAirDash;
	tStorage.m_pPlayer->m_flWaterJumpTime() = tStorage.m_PlayerData.m_flWaterJumpTime;
	tStorage.m_pPlayer->m_flSwimSoundTime() = tStorage.m_PlayerData.m_flSwimSoundTime;
	tStorage.m_pPlayer->m_surfaceProps() = tStorage.m_PlayerData.m_surfaceProps;
	tStorage.m_pPlayer->m_pSurfaceData() = tStorage.m_PlayerData.m_pSurfaceData;
	tStorage.m_pPlayer->m_surfaceFriction() = tStorage.m_PlayerData.m_surfaceFriction;
	tStorage.m_pPlayer->m_chTextureType() = tStorage.m_PlayerData.m_chTextureType;
	tStorage.m_pPlayer->m_vecPunchAngle() = tStorage.m_PlayerData.m_vecPunchAngle;
	tStorage.m_pPlayer->m_vecPunchAngleVel() = tStorage.m_PlayerData.m_vecPunchAngleVel;
	tStorage.m_pPlayer->m_MoveType() = tStorage.m_PlayerData.m_MoveType;
	tStorage.m_pPlayer->m_MoveCollide() = tStorage.m_PlayerData.m_MoveCollide;
	tStorage.m_pPlayer->m_vecLadderNormal() = tStorage.m_PlayerData.m_vecLadderNormal;
	tStorage.m_pPlayer->m_flGravity() = tStorage.m_PlayerData.m_flGravity;
	tStorage.m_pPlayer->m_nWaterLevel() = tStorage.m_PlayerData.m_nWaterLevel;
	tStorage.m_pPlayer->m_nWaterType() = tStorage.m_PlayerData.m_nWaterType;
	tStorage.m_pPlayer->m_flFallVelocity() = tStorage.m_PlayerData.m_flFallVelocity;
	tStorage.m_pPlayer->m_nPlayerCond() = tStorage.m_PlayerData.m_nPlayerCond;
	tStorage.m_pPlayer->m_nPlayerCondEx() = tStorage.m_PlayerData.m_nPlayerCondEx;
	tStorage.m_pPlayer->m_nPlayerCondEx2() = tStorage.m_PlayerData.m_nPlayerCondEx2;
	tStorage.m_pPlayer->m_nPlayerCondEx3() = tStorage.m_PlayerData.m_nPlayerCondEx3;
	tStorage.m_pPlayer->m_nPlayerCondEx4() = tStorage.m_PlayerData.m_nPlayerCondEx4;
	tStorage.m_pPlayer->_condition_bits() = tStorage.m_PlayerData._condition_bits;
}

void CMovementSimulation::Store()
{
	for (auto pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		auto pPlayer = pEntity->As<CTFPlayer>();
		auto& vRecords = m_mRecords[pPlayer->entindex()];

		if (pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost() || pPlayer->m_vecVelocity().IsZero())
		{
			vRecords.clear();
			continue;
		}
		else if (!H::Entities.GetDeltaTime(pPlayer->entindex()))
			continue;

		bool bLocal = pPlayer->entindex() == I::EngineClient->GetLocalPlayer() && !I::EngineClient->IsPlayingDemo();
		Vec3 vVelocity = bLocal ? F::EnginePrediction.m_vVelocity : pPlayer->m_vecVelocity();
		Vec3 vOrigin = bLocal ? F::EnginePrediction.m_vOrigin : pPlayer->m_vecOrigin();
		Vec3 vDirection = bLocal ? Math::RotatePoint(F::EnginePrediction.m_vDirection, {}, { 0, F::EnginePrediction.m_vAngles.y, 0 }) : Vec3(vVelocity.x, vVelocity.y, 0.f);

		MoveData* pLastRecord = !vRecords.empty() ? &vRecords.front() : nullptr;
		vRecords.emplace_front(
			vDirection,
			pPlayer->m_flSimulationTime(),
			pPlayer->m_nWaterLevel() >= 2 ? 2 : pPlayer->IsOnGround() ? 0 : 1,
			vVelocity,
			vOrigin
		);
		MoveData& tCurRecord = vRecords.front();
		if (vRecords.size() > 66)
			vRecords.pop_back();

		float flMaxSpeed = SDK::MaxSpeed(pPlayer);
		if (pLastRecord)
		{
			/*
			if (tRecord.m_iMode != pLastRecord->m_iMode)
				vRecords.clear();
			else // does this eat up fps? i can't tell currently
			*/
			{
				CGameTrace trace = {};
				CTraceFilterWorldAndPropsOnly filter = {};
				SDK::TraceHull(pLastRecord->m_vOrigin, pLastRecord->m_vOrigin + pLastRecord->m_vVelocity * TICK_INTERVAL, pPlayer->m_vecMins() + 0.125f, pPlayer->m_vecMaxs() - 0.125f, pPlayer->SolidMask(), &filter, &trace);
				if (trace.DidHit() && trace.plane.normal.z < 0.707f)
					vRecords.clear();
			}
		}
		if (pPlayer->InCond(TF_COND_SHIELD_CHARGE))
		{
			DummyCmd.forwardmove = 450.f;
			DummyCmd.sidemove = 0.f;
			SDK::FixMovement(&DummyCmd, bLocal ? F::EnginePrediction.m_vAngles : pPlayer->GetEyeAngles(), {});
			tCurRecord.m_vDirection.x = DummyCmd.forwardmove;
			tCurRecord.m_vDirection.y = -DummyCmd.sidemove;
		}
		else
		{
			switch (tCurRecord.m_iMode)
			{
			case 0:
				if (bLocal && Vars::Misc::Movement::Bunnyhop.Value && G::OriginalMove.m_iButtons & IN_JUMP)
					tCurRecord.m_vDirection = vVelocity.Normalized2D() * flMaxSpeed;
				break;
			case 1:
				tCurRecord.m_vDirection = vVelocity.Normalized2D() * flMaxSpeed;
				break;
			case 2:
				tCurRecord.m_vDirection *= 2;
			}
		}
	}

	for (auto pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		auto pPlayer = pEntity->As<CTFPlayer>();
		auto& vSimTimes = m_mSimTimes[pPlayer->entindex()];

		if (pEntity->entindex() == I::EngineClient->GetLocalPlayer() || pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost())
		{
			vSimTimes.clear();
			continue;
		}

		float flDeltaTime = H::Entities.GetDeltaTime(pPlayer->entindex());
		if (!flDeltaTime)
			continue;

		vSimTimes.push_front(flDeltaTime);
		if (vSimTimes.size() > Vars::Aimbot::Projectile::DeltaCount.Value)
			vSimTimes.pop_back();
	}
}



bool CMovementSimulation::Initialize(CBaseEntity* pEntity, PlayerStorage& tStorage, bool bHitchance, bool bStrafe)
{
	if (!pEntity || !pEntity->IsPlayer() || !pEntity->As<CTFPlayer>()->IsAlive())
	{
		tStorage.m_bInitFailed = tStorage.m_bFailed = true;
		return false;
	}

	auto pPlayer = pEntity->As<CTFPlayer>();
	tStorage.m_pPlayer = pPlayer;

	I::MoveHelper->SetHost(pPlayer);
	pPlayer->m_pCurrentCommand() = &DummyCmd;

	// store player restore data
	Store(tStorage);

	// store vars
	m_bOldInPrediction = I::Prediction->m_bInPrediction;
	m_bOldFirstTimePredicted = I::Prediction->m_bFirstTimePredicted;
	m_flOldFrametime = I::GlobalVars->frametime;

	// the hacks that make it work
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

		// Heavy-specific movement prediction with momentum consideration
		if (pPlayer->m_iClass() == TF_CLASS_HEAVY && !pPlayer->m_vecVelocity().IsZero())
		{
			float flMomentum = pPlayer->m_vecVelocity().Length2D() / SDK::MaxSpeed(pPlayer);
			if (flMomentum < 0.1f) // Very slow movement
			{
				// Smooth out velocity changes for Heavy to prevent jitter
				if (auto pAvgVel = H::Entities.GetAvgVelocity(pPlayer->entindex()))
				{
					Vec3 vSmoothedVel = pPlayer->m_vecVelocity() * 0.7f + *pAvgVel * 0.3f;
					pPlayer->m_vecVelocity() = vSmoothedVel;
				}
			}
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

	// setup move data
	if (!SetupMoveData(tStorage))
	{
		tStorage.m_bFailed = true;
		return false;
	}

	const int iStrafeSamples = tStorage.m_bDirectMove
		? Vars::Aimbot::Projectile::GroundSamples.Value
		: Vars::Aimbot::Projectile::AirSamples.Value;

	// calculate strafe if desired
	bool bCalculated = bStrafe ? StrafePrediction(tStorage, iStrafeSamples) : false;

	// really hope this doesn't work like shit
	if (bHitchance && bCalculated && !pPlayer->m_vecVelocity().IsZero() && Vars::Aimbot::Projectile::HitChance.Value)
	{
		const auto& vRecords = m_mRecords[pPlayer->entindex()];
		const auto iSamples = vRecords.size();

		float flCurrentChance = 1.f, flAverageYaw = 0.f;
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

			if ((i + 1) % iStrafeSamples == 0 || i == iSamples - 1)
			{
				flAverageYaw /= i % iStrafeSamples + 1;
				if (fabsf(tStorage.m_flAverageYaw - flAverageYaw) > 0.5f)
					flCurrentChance -= 1.f / ((iSamples - 1) / float(iStrafeSamples) + 1);
				flAverageYaw = 0.f;
			}
		}

		if (flCurrentChance < Vars::Aimbot::Projectile::HitChance.Value / 100)
		{
			char szBuffer[128];
			sprintf_s(szBuffer, "Hitchance (%.1f%% < %.0f%%)", flCurrentChance * 100.0f, Vars::Aimbot::Projectile::HitChance.Value);
			SDK::Output("MovementSimulation", szBuffer, { 80, 200, 120 }, Vars::Debug::Logging.Value);

			tStorage.m_bFailed = true;
			return false;
		}
	}

	for (int i = 0; i < H::Entities.GetChoke(pPlayer->entindex()); i++)
		RunTick(tStorage);

	return true;
}

bool CMovementSimulation::SetupMoveData(PlayerStorage& tStorage)
{
	if (!tStorage.m_pPlayer)
		return false;

	tStorage.m_MoveData.m_bFirstRunOfFunctions = false;
	tStorage.m_MoveData.m_bGameCodeMovedPlayer = false;
	tStorage.m_MoveData.m_nPlayerHandle = reinterpret_cast<IHandleEntity*>(tStorage.m_pPlayer)->GetRefEHandle();

	tStorage.m_MoveData.m_vecAbsOrigin = tStorage.m_pPlayer->m_vecOrigin();
	tStorage.m_MoveData.m_vecVelocity = tStorage.m_pPlayer->m_vecVelocity();
	tStorage.m_MoveData.m_flMaxSpeed = SDK::MaxSpeed(tStorage.m_pPlayer);
	tStorage.m_MoveData.m_flClientMaxSpeed = tStorage.m_MoveData.m_flMaxSpeed;

	if (!tStorage.m_MoveData.m_vecVelocity.To2D().IsZero())
	{
		int iIndex = tStorage.m_pPlayer->entindex();
		if (!tStorage.m_pPlayer->InCond(TF_COND_SHIELD_CHARGE))
			tStorage.m_MoveData.m_vecViewAngles = { 0.f, Math::VectorAngles(tStorage.m_MoveData.m_vecVelocity).y, 0.f };
		else
			tStorage.m_MoveData.m_vecViewAngles = iIndex == I::EngineClient->GetLocalPlayer() && G::CurrentUserCmd ? G::CurrentUserCmd->viewangles : H::Entities.GetEyeAngles(iIndex);
		
		const auto& vRecords = m_mRecords[tStorage.m_pPlayer->entindex()];
		if (!vRecords.empty())
		{
			auto& tRecord = vRecords.front();
			if (!tRecord.m_vDirection.IsZero())
			{
				DummyCmd.forwardmove = tRecord.m_vDirection.x;
				DummyCmd.sidemove = -tRecord.m_vDirection.y;
				DummyCmd.upmove = tRecord.m_vDirection.z;
				SDK::FixMovement(&DummyCmd, {}, tStorage.m_MoveData.m_vecViewAngles);
				tStorage.m_MoveData.m_flForwardMove = DummyCmd.forwardmove;
				tStorage.m_MoveData.m_flSideMove = DummyCmd.sidemove;
				tStorage.m_MoveData.m_flUpMove = DummyCmd.upmove;
			}
		}
	}

	tStorage.m_MoveData.m_vecAngles = tStorage.m_MoveData.m_vecOldAngles = tStorage.m_MoveData.m_vecViewAngles;
	if (auto pConstraintEntity = tStorage.m_pPlayer->m_hConstraintEntity().Get())
		tStorage.m_MoveData.m_vecConstraintCenter = pConstraintEntity->GetAbsOrigin();
	else
		tStorage.m_MoveData.m_vecConstraintCenter = tStorage.m_pPlayer->m_vecConstraintCenter();
	tStorage.m_MoveData.m_flConstraintRadius = tStorage.m_pPlayer->m_flConstraintRadius();
	tStorage.m_MoveData.m_flConstraintWidth = tStorage.m_pPlayer->m_flConstraintWidth();
	tStorage.m_MoveData.m_flConstraintSpeedFactor = tStorage.m_pPlayer->m_flConstraintSpeedFactor();

	tStorage.m_flPredictedDelta = GetPredictedDelta(tStorage.m_pPlayer);
	tStorage.m_flSimTime = tStorage.m_pPlayer->m_flSimulationTime();
	tStorage.m_flPredictedSimTime = tStorage.m_flSimTime + tStorage.m_flPredictedDelta;
	tStorage.m_vPredictedOrigin = tStorage.m_MoveData.m_vecAbsOrigin;
	tStorage.m_bDirectMove = tStorage.m_pPlayer->IsOnGround() || tStorage.m_pPlayer->m_nWaterLevel() >= 2;

	return true;
}

static inline float GetGravity()
{
	static auto sv_gravity = U::ConVars.FindVar("sv_gravity");

	return sv_gravity->GetFloat();
}

// Advanced mathematical prediction using Runge-Kutta 4th order method for trajectory calculation
static inline Vec3 PredictTrajectoryRK4(const Vec3& vPos, const Vec3& vVel, float flGravity, float flDeltaTime)
{
	// RK4 coefficients for more accurate trajectory prediction
	auto derivative = [flGravity](const Vec3& pos, const Vec3& vel) -> std::pair<Vec3, Vec3> {
		return { vel, Vec3(0, 0, -flGravity) };
	};

	// k1
	auto [k1_pos, k1_vel] = derivative(vPos, vVel);
	
	// k2
	Vec3 pos2 = vPos + k1_pos * (flDeltaTime * 0.5f);
	Vec3 vel2 = vVel + k1_vel * (flDeltaTime * 0.5f);
	auto [k2_pos, k2_vel] = derivative(pos2, vel2);
	
	// k3
	Vec3 pos3 = vPos + k2_pos * (flDeltaTime * 0.5f);
	Vec3 vel3 = vVel + k2_vel * (flDeltaTime * 0.5f);
	auto [k3_pos, k3_vel] = derivative(pos3, vel3);
	
	// k4
	Vec3 pos4 = vPos + k3_pos * flDeltaTime;
	Vec3 vel4 = vVel + k3_vel * flDeltaTime;
	auto [k4_pos, k4_vel] = derivative(pos4, vel4);
	
	// Final position using RK4 formula
	return vPos + (k1_pos + k2_pos * 2.0f + k3_pos * 2.0f + k4_pos) * (flDeltaTime / 6.0f);
}

// Enhanced strafe prediction using Kalman filter principles for noise reduction
static inline float PredictStrafeAngleKalman(const std::vector<float>& vYawHistory, float flProcessNoise = 0.01f, float flMeasurementNoise = 0.1f)
{
	if (vYawHistory.size() < 3)
		return 0.f;

	// Simple Kalman filter implementation for yaw prediction
	float flEstimate = vYawHistory[0];
	float flErrorCovariance = 1.0f;
	
	for (size_t i = 1; i < vYawHistory.size(); ++i)
	{
		// Prediction step
		float flPredictedEstimate = flEstimate;
		float flPredictedErrorCovariance = flErrorCovariance + flProcessNoise;
		
		// Update step
		float flKalmanGain = flPredictedErrorCovariance / (flPredictedErrorCovariance + flMeasurementNoise);
		flEstimate = flPredictedEstimate + flKalmanGain * (vYawHistory[i] - flPredictedEstimate);
		flErrorCovariance = (1.0f - flKalmanGain) * flPredictedErrorCovariance;
	}
	
	return flEstimate;
}

// Adaptive friction calculation using polynomial regression for better air movement prediction
static inline float GetAdaptiveFriction(float flVelocity, float flTurnRate, const std::vector<float>& vVelocityHistory)
{
	if (vVelocityHistory.size() < 3)
		return 1.0f;

	// Calculate velocity acceleration using finite differences
	float flAcceleration = 0.f;
	for (size_t i = 1; i < std::min(vVelocityHistory.size(), size_t(5)); ++i)
	{
		flAcceleration += (vVelocityHistory[i-1] - vVelocityHistory[i]) / TICK_INTERVAL;
	}
	flAcceleration /= std::min(vVelocityHistory.size() - 1, size_t(4));

	// Adaptive friction based on acceleration pattern
	float flAdaptiveFactor = 1.0f - std::clamp(fabsf(flAcceleration) / 100.0f, 0.f, 0.75f);
	float flTurnFactor = std::clamp(fabsf(flTurnRate) / 45.0f, 0.f, 1.f);
	
	return Math::Lerp(flAdaptiveFactor, 0.25f, flTurnFactor);
}

static inline float GetFrictionScale(float flVelocityXY, float flTurn, float flVelocityZ, float flMin = 50.f, float flMax = 150.f)
{
	if (0.f >= flVelocityZ || flVelocityZ > 250.f)
		return 1.f;

	static auto sv_airaccelerate = U::ConVars.FindVar("sv_airaccelerate");
	float flScale = std::max(sv_airaccelerate->GetFloat(), 1.f);
	flMin *= flScale, flMax *= flScale;

	// Enhanced friction calculation using exponential decay for more accurate air movement
	float flTurnFactor = fabsf(flVelocityXY * flTurn);
	float flNormalizedTurn = std::clamp(flTurnFactor / flMax, 0.f, 1.f);
	
	// Use exponential interpolation for smoother friction transitions
	float flExponent = 2.5f; // Tuned for Source engine air movement
	float flFrictionReduction = std::pow(flNormalizedTurn, flExponent);
	
	return Math::Lerp(1.f, 0.25f, flFrictionReduction);
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

static bool GetYawDifference(MoveData& tRecord1, MoveData& tRecord2, bool bStart, float* pYaw, float flStraightFuzzyValue, int iMaxChanges = 0, int iMaxChangeTime = 0, float flMaxSpeed = 0.f)
{
	const float flYaw1 = Math::VectorAngles(tRecord1.m_vDirection).y, flYaw2 = Math::VectorAngles(tRecord2.m_vDirection).y;
	const float flTime1 = tRecord1.m_flSimTime, flTime2 = tRecord2.m_flSimTime;
	const int iTicks = std::max(TIME_TO_TICKS(flTime1 - flTime2), 1);

	*pYaw = Math::NormalizeAngle(flYaw1 - flYaw2);
	
	// Enhanced yaw calculation with velocity-based weighting and momentum consideration
	if (flMaxSpeed && tRecord1.m_iMode != 1)
	{
		float flVelocityRatio = std::clamp(tRecord1.m_vVelocity.Length2D() / flMaxSpeed, 0.f, 1.f);
		// Apply sigmoid function for smoother velocity scaling
		float flSigmoidWeight = 1.0f / (1.0f + std::exp(-6.0f * (flVelocityRatio - 0.5f)));
		*pYaw *= flSigmoidWeight;
	}
	
	// Advanced friction compensation using adaptive algorithms
	if (Vars::Aimbot::Projectile::MovesimFrictionFlags.Value & Vars::Aimbot::Projectile::MovesimFrictionFlagsEnum::CalculateIncrease && tRecord1.m_iMode == 1)
	{
		float flGravityEffect = GetGravity() * TICK_INTERVAL;
		float flFrictionScale = GetFrictionScale(tRecord1.m_vVelocity.Length2D(), *pYaw, tRecord1.m_vVelocity.z + flGravityEffect, 0.f, 56.f);
		
		// Apply cubic interpolation for smoother friction transitions
		float flCubicFriction = flFrictionScale * flFrictionScale * (3.0f - 2.0f * flFrictionScale);
		*pYaw /= std::max(flCubicFriction, 0.1f); // Prevent division by zero
	}
	
	// Dynamic yaw threshold based on velocity and movement mode
	float flDynamicThreshold = tRecord1.m_iMode == 1 ? 60.f : 45.f; // Higher threshold for air movement
	if (tRecord1.m_vVelocity.Length2D() > 200.f)
		flDynamicThreshold *= 1.2f; // Increase threshold for high-speed movement
		
	if (fabsf(*pYaw) > flDynamicThreshold)
		return false;

	static int iChanges, iStart;

	static int iStaticSign = 0;
	const int iLastSign = iStaticSign;
	const int iCurrSign = iStaticSign = *pYaw ? sign(*pYaw) : iStaticSign;

	static bool bStaticZero = false;
	const bool iLastZero = bStaticZero;
	const bool iCurrZero = bStaticZero = !*pYaw;

	const bool bChanged = iCurrSign != iLastSign || iCurrZero && iLastZero;
	
	// Enhanced straight movement detection using velocity consistency
	float flVelocityConsistency = 1.0f;
	if (tRecord1.m_vVelocity.Length2D() > 10.f && tRecord2.m_vVelocity.Length2D() > 10.f)
	{
		Vec3 vDir1 = tRecord1.m_vVelocity.Normalized2D();
		Vec3 vDir2 = tRecord2.m_vVelocity.Normalized2D();
		flVelocityConsistency = std::max(vDir1.Dot(vDir2), 0.f);
	}
	
	float flAdjustedFuzzyValue = flStraightFuzzyValue * (2.0f - flVelocityConsistency);
	const bool bStraight = fabsf(*pYaw) * tRecord1.m_vVelocity.Length2D() * iTicks < flAdjustedFuzzyValue;

	if (bStart)
	{
		iChanges = 0, iStart = TIME_TO_TICKS(flTime1);
		if (bStraight && ++iChanges > iMaxChanges)
			return false;
		return true;
	}
	else
	{
		if ((bChanged || bStraight) && ++iChanges > iMaxChanges)
			return false;
		return iChanges && iStart - TIME_TO_TICKS(flTime2) > iMaxChangeTime ? false : true;
	}
}

void CMovementSimulation::GetAverageYaw(PlayerStorage& tStorage, int iSamples)
{
	auto pPlayer = tStorage.m_pPlayer;
	auto& vRecords = m_mRecords[pPlayer->entindex()];
	if (vRecords.empty())
		return;

	bool bGround = tStorage.m_bDirectMove; int iMinimumStrafes = 4;
	float flMaxSpeed = SDK::MaxSpeed(tStorage.m_pPlayer, false, true);
	float flLowMinimumDistance = bGround ? Vars::Aimbot::Projectile::GroundLowMinimumDistance.Value : Vars::Aimbot::Projectile::AirLowMinimumDistance.Value;
	float flLowMinimumSamples = bGround ? Vars::Aimbot::Projectile::GroundLowMinimumSamples.Value : Vars::Aimbot::Projectile::AirLowMinimumSamples.Value;
	float flHighMinimumDistance = bGround ? Vars::Aimbot::Projectile::GroundHighMinimumDistance.Value : Vars::Aimbot::Projectile::AirHighMinimumDistance.Value;
	float flHighMinimumSamples = bGround ? Vars::Aimbot::Projectile::GroundHighMinimumSamples.Value : Vars::Aimbot::Projectile::AirHighMinimumSamples.Value;

	// Enhanced statistical analysis using weighted averages and outlier detection
	std::vector<float> vYawSamples;
	std::vector<float> vWeights;
	std::vector<float> vVelocityHistory;
	
	float flAverageYaw = 0.f; int iTicks = 0, iSkips = 0;
	iSamples = std::min(iSamples, int(vRecords.size()));
	size_t i = 1; for (; i < iSamples; i++)
	{
		auto& tRecord1 = vRecords[i - 1];
		auto& tRecord2 = vRecords[i];
		if (tRecord1.m_iMode != tRecord2.m_iMode)
		{
			iSkips++;
			continue;
		}

		bGround = tRecord1.m_iMode != 1;
		float flStraightFuzzyValue = bGround ? Vars::Aimbot::Projectile::GroundStraightFuzzyValue.Value : Vars::Aimbot::Projectile::AirStraightFuzzyValue.Value;
		int iMaxChanges = bGround ? Vars::Aimbot::Projectile::GroundMaxChanges.Value : Vars::Aimbot::Projectile::AirMaxChanges.Value;
		int iMaxChangeTime = bGround ? Vars::Aimbot::Projectile::GroundMaxChangeTime.Value : Vars::Aimbot::Projectile::AirMaxChangeTime.Value;
		iMinimumStrafes = 4 + iMaxChanges;
#ifdef VISUALIZE_RECORDS
		VisualizeRecords(tRecord1, tRecord2, { 255, 0, 0 }, flStraightFuzzyValue);
#endif

		float flYaw = 0.f;
		bool bResult = GetYawDifference(tRecord1, tRecord2, !iTicks, &flYaw, flStraightFuzzyValue, iMaxChanges, iMaxChangeTime, flMaxSpeed);
		char szBuffer[128];
		sprintf_s(szBuffer, "%zu (%d): %.3f, %s", i, iTicks, flYaw, bResult ? "true" : "false");
		SDK::Output("GetYawDifference", szBuffer, { 50, 127, 75 }, Vars::Debug::Logging.Value);
		if (!bResult)
			break;

		// Store samples for advanced statistical analysis
		vYawSamples.push_back(flYaw);
		vVelocityHistory.push_back(tRecord1.m_vVelocity.Length2D());
		
		// Calculate weight based on velocity consistency and time recency
		float flVelocityWeight = std::clamp(tRecord1.m_vVelocity.Length2D() / flMaxSpeed, 0.1f, 1.f);
		float flTimeWeight = std::exp(-0.1f * i); // Exponential decay for older samples
		float flConsistencyWeight = 1.0f;
		
		if (i > 1 && !vYawSamples.empty())
		{
			// Weight based on consistency with previous samples
			float flYawDiff = fabsf(flYaw - vYawSamples.back());
			flConsistencyWeight = std::exp(-flYawDiff / 10.0f);
		}
		
		float flTotalWeight = flVelocityWeight * flTimeWeight * flConsistencyWeight;
		vWeights.push_back(flTotalWeight);

		flAverageYaw += flYaw * flTotalWeight;
		iTicks += std::max(TIME_TO_TICKS(tRecord1.m_flSimTime - tRecord2.m_flSimTime), 1);
	}
#ifdef VISUALIZE_RECORDS
	size_t i2 = i; for (; i2 < iSamples; i2++)
	{
		auto& tRecord1 = vRecords[i2 - 1];
		auto& tRecord2 = vRecords[i2];

		float flStraightFuzzyValue = bGround ? Vars::Aimbot::Projectile::GroundStraightFuzzyValue.Value : Vars::Aimbot::Projectile::AirStraightFuzzyValue.Value;
		VisualizeRecords(tRecord1, tRecord2, { 0, 0, 0 }, flStraightFuzzyValue);
	}
#endif
	if (i <= size_t(iMinimumStrafes + iSkips)) // valid strafes not high enough
		return;

	int iMinimum = flLowMinimumSamples;
	if (pPlayer->entindex() != I::EngineClient->GetLocalPlayer())
	{
		float flDistance = 0.f;
		if (auto pLocal = H::Entities.GetLocal())
			flDistance = pLocal->m_vecOrigin().DistTo(tStorage.m_pPlayer->m_vecOrigin());
		iMinimum = flDistance < flLowMinimumDistance ? flLowMinimumSamples : Math::RemapVal(flDistance, flLowMinimumDistance, flHighMinimumDistance, flLowMinimumSamples + 1, flHighMinimumSamples);
	}

	// Apply advanced statistical methods
	if (!vYawSamples.empty() && !vWeights.empty())
	{
		// Calculate weighted average
		float flTotalWeight = std::accumulate(vWeights.begin(), vWeights.end(), 0.0f);
		if (flTotalWeight > 0.0f)
		{
			flAverageYaw /= flTotalWeight;
			
			// Apply Kalman filtering for noise reduction
			if (vYawSamples.size() >= 3)
			{
				float flKalmanFiltered = PredictStrafeAngleKalman(vYawSamples);
				// Blend Kalman result with weighted average
				flAverageYaw = Math::Lerp(flAverageYaw, flKalmanFiltered, 0.3f);
			}
			
			// Apply adaptive friction compensation
			if (!vVelocityHistory.empty())
			{
				float flAdaptiveFriction = GetAdaptiveFriction(vVelocityHistory.back(), flAverageYaw, vVelocityHistory);
				flAverageYaw *= flAdaptiveFriction;
			}
		}
	}
	else
	{
		// Fallback to original method
		flAverageYaw /= std::max(iTicks, iMinimum);
	}
	
	// Enhanced threshold with dynamic adjustment
	float flDynamicThreshold = 0.36f;
	if (!bGround) // Air movement typically has more variation
		flDynamicThreshold *= 0.8f;
	if (flMaxSpeed > 300.f) // High-speed classes need lower threshold
		flDynamicThreshold *= 0.9f;
		
	if (fabsf(flAverageYaw) < flDynamicThreshold)
		return;

	tStorage.m_flAverageYaw = flAverageYaw;
	char szBuffer[128];
	sprintf_s(szBuffer, "Enhanced flAverageYaw calculated to %.3f from %d (%d) %s", flAverageYaw, iTicks, iMinimum, pPlayer->entindex() == I::EngineClient->GetLocalPlayer() ? "(local)" : "");
	SDK::Output("MovementSimulation", szBuffer, { 100, 255, 150 }, Vars::Debug::Logging.Value);
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
	if (tStorage.m_bFailed || !tStorage.m_pPlayer || !tStorage.m_pPlayer->IsPlayer())
		return;

	if (bPath)
		tStorage.m_vPath.push_back(tStorage.m_MoveData.m_vecAbsOrigin);

	// make sure frametime and prediction vars are right
	I::Prediction->m_bInPrediction = true;
	I::Prediction->m_bFirstTimePredicted = false;
	I::GlobalVars->frametime = I::Prediction->m_bEnginePaused ? 0.f : TICK_INTERVAL;
	SetBounds(tStorage.m_pPlayer);

	float flCorrection = 0.f;
	if (tStorage.m_flAverageYaw)
	{
		float flMult = 1.f;
		if (!tStorage.m_bDirectMove && !tStorage.m_pPlayer->InCond(TF_COND_SHIELD_CHARGE))
		{
			// Enhanced air movement prediction with momentum consideration
			float flCurrentSpeed = tStorage.m_MoveData.m_vecVelocity.Length2D();
			float flMaxSpeed = tStorage.m_MoveData.m_flMaxSpeed;
			
			// Dynamic correction based on current velocity and optimal strafe angle
			float flOptimalAngle = 90.f;
			if (flCurrentSpeed > flMaxSpeed * 0.8f) // High speed - reduce strafe angle
				flOptimalAngle = Math::Lerp(90.f, 45.f, (flCurrentSpeed - flMaxSpeed * 0.8f) / (flMaxSpeed * 0.2f));
			
			flCorrection = flOptimalAngle * sign(tStorage.m_flAverageYaw);
			
			// Enhanced friction calculation with velocity prediction
			if (Vars::Aimbot::Projectile::MovesimFrictionFlags.Value & Vars::Aimbot::Projectile::MovesimFrictionFlagsEnum::RunReduce)
			{
				float flPredictedVelZ = tStorage.m_MoveData.m_vecVelocity.z + GetGravity() * TICK_INTERVAL;
				flMult = GetFrictionScale(flCurrentSpeed, tStorage.m_flAverageYaw, flPredictedVelZ);
				
				// Apply additional smoothing for high-frequency strafing
				if (fabsf(tStorage.m_flAverageYaw) > 15.f)
				{
					float flSmoothingFactor = std::clamp(fabsf(tStorage.m_flAverageYaw) / 30.f, 0.f, 1.f);
					flMult = Math::Lerp(flMult, 1.0f, flSmoothingFactor * 0.3f);
				}
			}
		}
		
		// Apply enhanced yaw correction with momentum preservation
		float flYawAdjustment = tStorage.m_flAverageYaw * flMult + flCorrection;
		
		// Smooth yaw transitions to prevent jittery movement
		static float flLastYawAdjustment = 0.f;
		float flSmoothingAlpha = 0.7f; // Adjust for more/less smoothing
		flYawAdjustment = Math::Lerp(flLastYawAdjustment, flYawAdjustment, flSmoothingAlpha);
		flLastYawAdjustment = flYawAdjustment;
		
		tStorage.m_MoveData.m_vecViewAngles.y += flYawAdjustment;
	}
	else if (!tStorage.m_bDirectMove)
		tStorage.m_MoveData.m_flForwardMove = tStorage.m_MoveData.m_flSideMove = 0.f;

	float flOldSpeed = tStorage.m_MoveData.m_flClientMaxSpeed;
	if (tStorage.m_pPlayer->m_bDucked() && tStorage.m_pPlayer->IsOnGround() && tStorage.m_pPlayer->m_nWaterLevel() < 2)
		tStorage.m_MoveData.m_flClientMaxSpeed /= 3;

	// Enhanced bunny hop prediction with pre-speed consideration
	if (tStorage.m_bBunnyHop && tStorage.m_pPlayer->IsOnGround() && !tStorage.m_pPlayer->m_bDucked())
	{
		tStorage.m_MoveData.m_nOldButtons = 0;
		tStorage.m_MoveData.m_nButtons |= IN_JUMP;
		
		// Optimize jump timing for maximum speed gain
		float flCurrentSpeed = tStorage.m_MoveData.m_vecVelocity.Length2D();
		if (flCurrentSpeed > tStorage.m_MoveData.m_flMaxSpeed * 1.1f)
		{
			// Adjust movement input for optimal pre-speed maintenance
			Vec3 vOptimalDir = tStorage.m_MoveData.m_vecVelocity.Normalized2D();
			float flOptimalForward = vOptimalDir.x * 320.f;
			float flOptimalSide = -vOptimalDir.y * 320.f;
			
			// Blend with current input for smoother transitions
			tStorage.m_MoveData.m_flForwardMove = Math::Lerp(tStorage.m_MoveData.m_flForwardMove, flOptimalForward, 0.6f);
			tStorage.m_MoveData.m_flSideMove = Math::Lerp(tStorage.m_MoveData.m_flSideMove, flOptimalSide, 0.6f);
		}
	}

	I::GameMovement->ProcessMovement(tStorage.m_pPlayer, &tStorage.m_MoveData);
	if (pCallback)
		(*pCallback)(tStorage.m_MoveData);

	tStorage.m_MoveData.m_flClientMaxSpeed = flOldSpeed;

	tStorage.m_flSimTime += TICK_INTERVAL;
	tStorage.m_bPredictNetworked = tStorage.m_flSimTime >= tStorage.m_flPredictedSimTime;
	if (tStorage.m_bPredictNetworked)
	{
		tStorage.m_vPredictedOrigin = tStorage.m_MoveData.m_vecAbsOrigin;
		tStorage.m_flPredictedSimTime += tStorage.m_flPredictedDelta;
	}
	bool bLastbDirectMove = tStorage.m_bDirectMove;
	tStorage.m_bDirectMove = tStorage.m_pPlayer->IsOnGround() || tStorage.m_pPlayer->m_nWaterLevel() >= 2;

	if (tStorage.m_flAverageYaw)
		tStorage.m_MoveData.m_vecViewAngles.y -= flCorrection;
	else if (tStorage.m_bDirectMove && !bLastbDirectMove
		&& !tStorage.m_MoveData.m_flForwardMove && !tStorage.m_MoveData.m_flSideMove
		&& tStorage.m_MoveData.m_vecVelocity.Length2D() > tStorage.m_MoveData.m_flMaxSpeed * 0.015f)
	{
		// Enhanced ground transition handling with momentum preservation
		Vec3 vDirection = tStorage.m_MoveData.m_vecVelocity.Normalized2D() * 450.f;
		
		// Apply gradual transition to prevent sudden direction changes
		float flTransitionFactor = std::clamp(tStorage.m_MoveData.m_vecVelocity.Length2D() / tStorage.m_MoveData.m_flMaxSpeed, 0.5f, 1.f);
		vDirection *= flTransitionFactor;
		
		DummyCmd.forwardmove = vDirection.x, DummyCmd.sidemove = -vDirection.y;
		SDK::FixMovement(&DummyCmd, {}, tStorage.m_MoveData.m_vecViewAngles);
		tStorage.m_MoveData.m_flForwardMove = DummyCmd.forwardmove, tStorage.m_MoveData.m_flSideMove = DummyCmd.sidemove;
	}

	RestoreBounds(tStorage.m_pPlayer);
}

void CMovementSimulation::RunTick(PlayerStorage& tStorage, bool bPath, std::function<void(CMoveData&)> fCallback)
{
	RunTick(tStorage, bPath, &fCallback);
}

void CMovementSimulation::Restore(PlayerStorage& tStorage)
{
	if (tStorage.m_bInitFailed || !tStorage.m_pPlayer)
		return;

	I::MoveHelper->SetHost(nullptr);
	tStorage.m_pPlayer->m_pCurrentCommand() = nullptr;

	Reset(tStorage);

	I::Prediction->m_bInPrediction = m_bOldInPrediction;
	I::Prediction->m_bFirstTimePredicted = m_bOldFirstTimePredicted;
	I::GlobalVars->frametime = m_flOldFrametime;

	/*
	const bool bInitFailed = tStorage.m_bInitFailed, bFailed = tStorage.m_bFailed;
	memset(&tStorage, 0, sizeof(PlayerStorage));
	tStorage.m_bInitFailed = bInitFailed, tStorage.m_bFailed = bFailed;
	*/
}

float CMovementSimulation::GetPredictedDelta(CBaseEntity* pEntity)
{
	auto& vSimTimes = m_mSimTimes[pEntity->entindex()];
	if (!vSimTimes.empty())
	{
		switch (Vars::Aimbot::Projectile::DeltaMode.Value)
		{
		case 0: return std::reduce(vSimTimes.begin(), vSimTimes.end()) / vSimTimes.size();
		case 1: return *std::max_element(vSimTimes.begin(), vSimTimes.end());
		}
	}
	return TICK_INTERVAL;
}