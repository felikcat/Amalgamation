// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "AimbotMelee.h"

#include "../Aimbot.h"
#include "../../Simulation/MovementSimulation/MovementSimulation.h"
#include "../../Ticks/Ticks.h"
#include "../../Visuals/Visuals.h"

std::vector<Target_t> CAimbotMelee::GetTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon)
{
	std::vector<Target_t> vTargets;
	vTargets.reserve(64); // Reserve space to avoid frequent reallocations

	const Vec3 vLocalPos = F::Ticks.GetShootPos();
	const Vec3 vLocalAngles = I::EngineClient->GetViewAngles();

	if (Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::Players)
	{
		const bool bDisciplinary = Vars::Aimbot::Melee::WhipTeam.Value && SDK::AttribHookValue(0, "speed_buff_ally", pWeapon) > 0;
		const auto& entityGroup = H::Entities.GetGroup(bDisciplinary ? EGroupType::PLAYERS_ALL : EGroupType::PLAYERS_ENEMIES);
		const int localTeam = pLocal->m_iTeamNum();
		
		for (const auto pEntity : entityGroup)
		{
			const bool bTeammate = pEntity->m_iTeamNum() == localTeam;
			if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
				continue;

			float flFOVTo;
			Vec3 vPos, vAngleTo;
			if (!F::AimbotGlobal.PlayerBoneInFOV(pEntity->As<CTFPlayer>(), vLocalPos, vLocalAngles, flFOVTo, vPos, vAngleTo))
				continue;

			const float flDistTo = vLocalPos.DistTo(vPos);
			vTargets.emplace_back(pEntity, TargetEnum::Player, vPos, vAngleTo, flFOVTo, flDistTo,
								  bTeammate ? 0 : F::AimbotGlobal.GetPriority(pEntity->entindex()));
		}
	}

	if (Vars::Aimbot::General::Target.Value)
	{
		const bool bWrench = pWeapon->GetWeaponID() == TF_WEAPON_WRENCH;
		const bool bDestroySapper = pWeapon->GetWeaponID() == TF_WEAPON_FIREAXE && SDK::AttribHookValue(0, "set_dmg_apply_to_sapper", pWeapon);
		const auto& buildingGroup = H::Entities.GetGroup(bWrench || bDestroySapper ? EGroupType::BUILDINGS_ALL : EGroupType::BUILDINGS_ENEMIES);
		const int localTeam = pLocal->m_iTeamNum();
		const float maxFOV = Vars::Aimbot::General::AimFOV.Value;

		for (const auto pEntity : buildingGroup)
		{
			if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
				continue;

			if (pEntity->m_iTeamNum() == localTeam &&
				((bWrench && !AimFriendlyBuilding(pEntity->As<CBaseObject>())) ||
				 (bDestroySapper && !pEntity->As<CBaseObject>()->m_bHasSapper())))
				continue;

			const Vec3 vPos = pEntity->GetCenter();
			const Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPos);
			const float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
			if (flFOVTo > maxFOV)
				continue;

			const float flDistTo = vLocalPos.DistTo(vPos);
			const int targetType = pEntity->IsSentrygun() ? TargetEnum::Sentry :
								   pEntity->IsDispenser() ? TargetEnum::Dispenser : TargetEnum::Teleporter;
			vTargets.emplace_back(pEntity, targetType, vPos, vAngleTo, flFOVTo, flDistTo);
		}
	}

	if (Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::NPCs)
	{
		const auto& npcGroup = H::Entities.GetGroup(EGroupType::WORLD_NPC);
		const float maxFOV = Vars::Aimbot::General::AimFOV.Value;

		for (const auto pEntity : npcGroup)
		{
			if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
				continue;

			const Vec3 vPos = pEntity->GetCenter();
			const Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPos);
			const float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
			if (flFOVTo > maxFOV)
				continue;

			const float flDistTo = vLocalPos.DistTo(vPos);
			vTargets.emplace_back(pEntity, TargetEnum::NPC, vPos, vAngleTo, flFOVTo, flDistTo);
		}
	}

	return vTargets;
}

bool CAimbotMelee::AimFriendlyBuilding(CBaseObject* pBuilding)
{
	if (!pBuilding->m_bMiniBuilding() && pBuilding->m_iUpgradeLevel() != 3 || pBuilding->m_iHealth() < pBuilding->m_iMaxHealth() || pBuilding->m_bHasSapper())
		return true;

	if (pBuilding->IsSentrygun())
	{
		int iShells, iMaxShells, iRockets, iMaxRockets; pBuilding->As<CObjectSentrygun>()->GetAmmoCount(iShells, iMaxShells, iRockets, iMaxRockets);
		if (iShells < iMaxShells || iRockets < iMaxRockets)
			return true;
	}

	return false;
}

std::vector<Target_t> CAimbotMelee::SortTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon)
{
	auto vTargets = GetTargets(pLocal, pWeapon);

	F::AimbotGlobal.SortTargets(vTargets, Vars::Aimbot::General::TargetSelectionEnum::Distance);
	
	const auto maxTargets = static_cast<std::size_t>(Vars::Aimbot::General::MaxTargets.Value);
	if (vTargets.size() > maxTargets) {
		vTargets.resize(maxTargets);
	}
	
	F::AimbotGlobal.SortPriority(vTargets);
	return vTargets;
}



int CAimbotMelee::GetSwingTime(CTFWeaponBase* pWeapon)
{
	if (pWeapon->GetWeaponID() == TF_WEAPON_KNIFE)
		return 0;
	return Vars::Aimbot::Melee::SwingTicks.Value;
}

void CAimbotMelee::SimulatePlayers(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, std::vector<Target_t> vTargets,
								   Vec3& vEyePos, std::unordered_map<int, std::deque<TickRecord>>& mRecordMap,
								   std::unordered_map<int, std::vector<Vec3>>& mPaths)
{
	// swing prediction / auto warp
	const int iSwingTicks = GetSwingTime(pWeapon);
	const int iMax = (iDoubletapTicks && Vars::Doubletap::AntiWarp.Value && pLocal->m_hGroundEntity())
		? std::max(iSwingTicks - Vars::Doubletap::TickLimit.Value - 1, 0)
		: std::max(iSwingTicks, iDoubletapTicks);

	if ((Vars::Aimbot::Melee::SwingPrediction.Value || iDoubletapTicks) && pWeapon->m_flSmackTime() < 0.f && iMax > 0)
	{
		PlayerStorage tStorage;
		std::unordered_map<int, PlayerStorage> mStorage;
		mStorage.reserve(vTargets.size()); // Reserve space for better performance

		F::MoveSim.Initialize(pLocal, tStorage, false, !iDoubletapTicks);
		for (const auto& tTarget : vTargets) {
			F::MoveSim.Initialize(tTarget.m_pEntity, mStorage[tTarget.m_pEntity->entindex()], false);
		}

		const int swingTicksMinusDoubletap = iSwingTicks - iDoubletapTicks;
		const int swingTime = GetSwingTime(pWeapon);
		const bool useSwingPredictLag = Vars::Aimbot::Melee::SwingPredictLag.Value;
		
		for (int i = 0; i < iMax; ++i) // intended for plocal to collide with targets
		{
			// Handle local player simulation
			if (pLocal->InCond(TF_COND_SHIELD_CHARGE) && (iMax - i) <= swingTime) {
				// demo charge fix for swing pred
				const float maxSpeed = SDK::MaxSpeed(pLocal, false, true);
				tStorage.m_MoveData.m_flMaxSpeed = tStorage.m_MoveData.m_flClientMaxSpeed = maxSpeed;
			}
			F::MoveSim.RunTick(tStorage);
			
			// Handle target simulation
			if (i < swingTicksMinusDoubletap)
			{
				for (const auto& tTarget : vTargets)
				{
					auto& targetStorage = mStorage[tTarget.m_pEntity->entindex()];

					F::MoveSim.RunTick(targetStorage);
					if (!targetStorage.m_bFailed) {
						const float simTime = (!useSwingPredictLag || targetStorage.m_bPredictNetworked)
							? tTarget.m_pEntity->m_flSimulationTime() + TICKS_TO_TIME(i + 1)
							: 0.f;
						const Vec3& origin = useSwingPredictLag
							? targetStorage.m_vPredictedOrigin
							: targetStorage.m_MoveData.m_vecAbsOrigin;
							
						mRecordMap[tTarget.m_pEntity->entindex()].emplace_front(
							simTime,
							origin,
							tTarget.m_pEntity->m_vecMins(),
							tTarget.m_pEntity->m_vecMaxs(),
							BoneMatrix{},
							false,
							origin
						);
					}
				}
			}
		}
		vEyePos = tStorage.m_MoveData.m_vecAbsOrigin + pLocal->m_vecViewOffset();

		if (Vars::Visuals::Simulation::SwingLines.Value && Vars::Visuals::Simulation::PlayerPath.Value)
		{
			const bool bAlwaysDraw = !Vars::Aimbot::General::AutoShoot.Value || Vars::Debug::Info.Value;
			if (!bAlwaysDraw)
			{
				mPaths[pLocal->entindex()] = tStorage.m_vPath;
				for (auto& tTarget : vTargets)
					mPaths[tTarget.m_pEntity->entindex()] = mStorage[tTarget.m_pEntity->entindex()].m_vPath;
			}
			else
			{
				G::LineStorage.clear();
				G::BoxStorage.clear();
				G::PathStorage.clear();

				if (Vars::Colors::PlayerPath.Value.a)
				{
					G::PathStorage.emplace_back(tStorage.m_vPath, I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
					for (auto& tTarget : vTargets)
						G::PathStorage.emplace_back(mStorage[tTarget.m_pEntity->entindex()].m_vPath, I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
				}
				if (Vars::Colors::PlayerPathClipped.Value.a)
				{
					G::PathStorage.emplace_back(tStorage.m_vPath, I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
					for (auto& tTarget : vTargets)
						G::PathStorage.emplace_back(mStorage[tTarget.m_pEntity->entindex()].m_vPath, I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
				}
			}
		}

		F::MoveSim.Restore(tStorage);
		for (auto& tTarget : vTargets)
			F::MoveSim.Restore(mStorage[tTarget.m_pEntity->entindex()]);
	}
}

// Helper function to normalize yaw angle to [-180, 180] range
static constexpr float NormalizeYaw(float yaw) noexcept {
	while (yaw > 180.0f) yaw -= 360.0f;
	while (yaw < -180.0f) yaw += 360.0f;
	return yaw;
}

// Simplified ping compensation helper - predicts target rotation
static Vec3 PredictTargetAngles(CBaseEntity* pTarget, float flPingSeconds) noexcept {
	if (!pTarget || flPingSeconds <= 0.0f) {
		return { 0.f, H::Entities.GetEyeAngles(pTarget->entindex()).y, 0.f };
	}

	// Get current and ping-compensated angles
	const Vec3 currentAngles = H::Entities.GetEyeAngles(pTarget->entindex());
	const Vec3 pingAngles = H::Entities.GetPingAngles(pTarget->entindex());
	
	// Use a more conservative prediction approach
	// Apply ping compensation with reduced intensity to avoid over-correction
	const float compensationFactor = std::clamp(flPingSeconds * 0.5f, 0.0f, 1.0f);
	const float predictedYaw = currentAngles.y + (pingAngles.y * compensationFactor);
	
	return { 0.f, NormalizeYaw(predictedYaw), 0.f };
}

// Simplified relative velocity compensation
static float CalculateRelativeVelocityScale(CTFPlayer* pLocal, CBaseEntity* pTarget, float flDist) noexcept {
	const Vec3 localVel = pLocal->m_vecVelocity();
	const float localSpeed = localVel.Length2D();
	
	// Conservative velocity-based scaling - only apply minor adjustments
	const float speedFactor = 1.0f + std::clamp(localSpeed / 400.0f, 0.0f, 0.25f); // Max 25% increase
	
	// Minimal distance scaling to avoid over-compensation
	const float distanceFactor = std::clamp(80.0f / std::max(flDist, 40.0f), 0.9f, 1.1f); // Max 10% adjustment
	
	return speedFactor * distanceFactor;
}

// Simplified dual-angle testing for high ping scenarios
static bool TestDualAngles(const std::function<bool(const Vec3&)>& testFunc,
						   CBaseEntity* pTarget, const Vec3& baseAngles,
						   float flPingSeconds) noexcept {
	// Test base angle first
	if (testFunc(baseAngles)) {
		return true;
	}

	// For high ping, test with conservative ping compensation
	if (flPingSeconds > 0.05f) { // Only for ping > 50ms
		const Vec3 pingAngles = H::Entities.GetPingAngles(pTarget->entindex());
		
		// Apply conservative compensation (max 50% of ping angle)
		const float compensationFactor = std::clamp(flPingSeconds, 0.0f, 0.5f);
		Vec3 compensatedAngles = baseAngles;
		compensatedAngles.y = NormalizeYaw(baseAngles.y + (pingAngles.y * compensationFactor));
		
		return testFunc(compensatedAngles);
	}
	
	return false;
}


// Improved backstab detection with simplified, more reliable compensation system
// Fixes issues with over-aggressive ping compensation that caused backstab failures
bool CAimbotMelee::CanBackstab(CBaseEntity* pTarget, CTFPlayer* pLocal, Vec3 vEyeAngles)
{
	if (!pLocal || !pTarget) {
		return false;
	}

	// Check for razorback protection
	if (Vars::Aimbot::Melee::IgnoreRazorback.Value) {
		CUtlVector<CBaseEntity*> itemList;
		const int iBackstabShield = SDK::AttribHookValue(0, "set_blockbackstab_once", pTarget, &itemList);
		if (iBackstabShield && itemList.Count() > 0) {
			CBaseEntity* pEntity = itemList.Element(0);
			if (pEntity && pEntity->ShouldDraw()) {
				return false;
			}
		}
	}

	// Calculate 2D vector from local player to target with improved precision
	const Vec3 vLocalOrigin = pLocal->m_vecOrigin();
	const Vec3 vTargetOrigin = pTarget->GetAbsOrigin();
	Vec3 vToTarget = (vTargetOrigin - vLocalOrigin).To2D();
	
	const float flDist = vToTarget.Normalize();
	if (flDist <= std::numeric_limits<float>::epsilon() * 10.0f) { // More precise epsilon check
		return false;
	}

	// Simplified ping compensation system
	const float flPing = I::EngineClient->GetNetChannelInfo()->GetLatency(FLOW_OUTGOING) * 1000.0f; // Convert to ms
	const float flPingSeconds = flPing / 1000.0f; // Convert back to seconds for calculations
	
	// Conservative tolerance calculations
	constexpr float flBaseTolerance = 0.0625f;     // Base network tolerance (1/16)
	constexpr float flMinDistance = 1.0f;          // Minimum distance for calculations
	constexpr float flMaxDistance = 128.0f;        // Maximum effective distance
	
	// Simplified ping scaling - much more conservative
	const float flPingScale = 1.0f + std::clamp(flPing / 200.0f, 0.0f, 1.0f); // Max 2x scaling for 200ms+ ping
	
	// Simplified velocity compensation
	const Vec3 localVelocity = pLocal->m_vecVelocity();
	const float localSpeed = localVelocity.Length2D();
	const float relativeVelocityScale = CalculateRelativeVelocityScale(pLocal, pTarget, flDist);
	
	// Conservative velocity scaling
	const float velocityScale = 1.0f + std::clamp(localSpeed / 400.0f, 0.0f, 0.3f); // Max 30% increase
	
	// Combined scaling with limits to prevent over-compensation
	const float combinedScale = std::clamp(flPingScale * velocityScale * relativeVelocityScale, 1.0f, 2.0f);
	
	// Calculate distance-based tolerance
	const float flClampedDistance = std::clamp(flDist, flMinDistance, flMaxDistance);
	const float flDistanceTolerance = (flBaseTolerance * combinedScale) / flClampedDistance;
	
	// TF2-specific backstab angle thresholds - keep original values for consistency
	constexpr float flBaseTargetViewDot = 0.0f;           // Target must be facing away (180° cone)
	constexpr float flBaseOwnerViewDot = 0.5f;            // ~60° cone behind target (cos(60°) = 0.5)
	constexpr float flBaseFacestabDot = -0.3f;            // Anti-facestab protection (~107° max)
	
	// Conservative network tolerance scaling
	const float flNetworkTolerance = 0.00306796f * std::clamp(combinedScale, 1.0f, 1.5f);
	
	// Apply tolerances with conservative adjustments
	const float flPosVsTargetViewMinDot = flBaseTargetViewDot + flNetworkTolerance + flDistanceTolerance;
	const float flPosVsOwnerViewMinDot = flBaseOwnerViewDot + (flDistanceTolerance * 0.5f); // Reduced impact
	const float flViewAnglesMinDot = flBaseFacestabDot + (flNetworkTolerance * 0.5f); // Reduced impact

	const auto TestDots = [&](const Vec3& vTargetAngles) -> bool {
		// Calculate forward vectors with improved precision
		Vec3 vOwnerForward;
		Math::AngleVectors(vEyeAngles, &vOwnerForward);
		
		Vec3 vTargetForward;
		Math::AngleVectors(vTargetAngles, &vTargetForward);
		
		// Normalize 2D vectors and validate
		const float ownerLength2D = vOwnerForward.Length2D();
		const float targetLength2D = vTargetForward.Length2D();
		
		if (ownerLength2D < 0.001f || targetLength2D < 0.001f) {
			return false; // Invalid vectors
		}
		
		vOwnerForward.x /= ownerLength2D;
		vOwnerForward.y /= ownerLength2D;
		vTargetForward.x /= targetLength2D;
		vTargetForward.y /= targetLength2D;

		// Calculate dot products for backstab validation
		const float flPosVsTargetViewDot = vToTarget.Dot(vTargetForward);   // Is attacker behind target?
		const float flPosVsOwnerViewDot = vToTarget.Dot(vOwnerForward);     // Is attacker facing target?
		const float flViewAnglesDot = vTargetForward.Dot(vOwnerForward);    // Anti-facestab check

		// Apply backstab validation with consistent thresholds
		const bool bBehindTarget = flPosVsTargetViewDot > flPosVsTargetViewMinDot;
		const bool bFacingTarget = flPosVsOwnerViewDot > flPosVsOwnerViewMinDot;
		const bool bNotFacestab = flViewAnglesDot > flViewAnglesMinDot;

		return bBehindTarget && bFacingTarget && bNotFacestab;
	};

	// Get current target angles (only yaw matters for backstab)
	Vec3 vTargetAngles{ 0.f, NormalizeYaw(H::Entities.GetEyeAngles(pTarget->entindex()).y), 0.f };
	
	const bool accountPing = Vars::Aimbot::Melee::BackstabAccountPing.Value;
	const bool doubleTest = Vars::Aimbot::Melee::BackstabDoubleTest.Value;
	
	if (!accountPing) {
		// Simple test without ping compensation
		return TestDots(vTargetAngles);
	}
	
	// Advanced ping compensation system
	const bool bHighPing = flPing > 100.0f;
	const bool bVeryHighPing = flPing > 150.0f;
	
	// For very high ping, use simplified dual-angle testing
	if (bVeryHighPing) {
		return TestDualAngles(TestDots, pTarget, vTargetAngles, flPingSeconds);
	}
	
	// Simplified high-ping compensation
	if (bHighPing) {
		// Test current angles first
		const bool bCurrentTest = TestDots(vTargetAngles);
		if (bCurrentTest) {
			return true;
		}
		
		// Test with conservative ping compensation
		const Vec3 vPingAngles = H::Entities.GetPingAngles(pTarget->entindex());
		Vec3 vPingCompensatedAngles = vTargetAngles;
		
		// Conservative ping multiplier - avoid over-compensation
		const float pingMultiplier = std::clamp(1.0f + (flPing - 100.0f) / 300.0f, 1.0f, 1.5f); // Max 1.5x
		vPingCompensatedAngles.y = NormalizeYaw(vTargetAngles.y + (vPingAngles.y * pingMultiplier));
		
		return TestDots(vPingCompensatedAngles);
	}
	
	// Standard ping compensation for moderate ping
	if (doubleTest) {
		// Test both current and ping-compensated angles (AND logic)
		const bool bCurrentTest = TestDots(vTargetAngles);
		if (!bCurrentTest) {
			return false; // Failed current test
		}
		
		const Vec3 vPingAngles = H::Entities.GetPingAngles(pTarget->entindex());
		Vec3 vPingCompensatedAngles = vTargetAngles;
		vPingCompensatedAngles.y = NormalizeYaw(vTargetAngles.y + vPingAngles.y);
		
		return TestDots(vPingCompensatedAngles);
	}
	
	// Single test with ping compensation (OR logic)
	const bool bCurrentTest = TestDots(vTargetAngles);
	if (bCurrentTest) {
		return true;
	}
	
	const Vec3 vPingAngles = H::Entities.GetPingAngles(pTarget->entindex());
	Vec3 vPingCompensatedAngles = vTargetAngles;
	vPingCompensatedAngles.y = NormalizeYaw(vTargetAngles.y + vPingAngles.y);
	
	return TestDots(vPingCompensatedAngles);
}

int CAimbotMelee::CanHit(Target_t& tTarget, CTFPlayer* pLocal, CTFWeaponBase* pWeapon, Vec3 vEyePos, std::deque<TickRecord>& vSimRecords)
{
	// Early return for unsimulated targets
	if (Vars::Aimbot::General::Ignore.Value & Vars::Aimbot::General::IgnoreEnum::Unsimulated &&
		H::Entities.GetChoke(tTarget.m_pEntity->entindex()) > Vars::Aimbot::General::TickTolerance.Value) {
		return false;
	}

	// Calculate weapon parameters
	const float baseHull = SDK::AttribHookValue(18, "melee_bounds_multiplier", pWeapon);
	const float baseRange = SDK::AttribHookValue(pWeapon->GetSwingRange(pLocal), "melee_range_multiplier", pWeapon);
	const float modelScale = pLocal->m_flModelScale();
	
	const float flHull = (modelScale > 1.0f) ? baseHull * modelScale : baseHull;
	const float flRange = (modelScale > 1.0f) ? baseRange * modelScale : baseRange;
	
	const Vec3 vSwingMins{ -flHull, -flHull, -flHull };
	const Vec3 vSwingMaxs{ flHull, flHull, flHull };

	std::vector<TickRecord*> vRecords;
	vRecords.reserve(32); // Reserve space for better performance
	
	if (F::Backtrack.GetRecords(tTarget.m_pEntity, vRecords)) {
		if (!vRecords.empty()) {
			// Add simulation records
			vRecords.reserve(vRecords.size() + vSimRecords.size());
			for (auto& tRecord : vSimRecords) {
				vRecords.push_back(&tRecord);
			}
			vRecords = F::Backtrack.GetValidRecords(vRecords, pLocal, true, -TICKS_TO_TIME(vSimRecords.size()));
		}
		if (vRecords.empty()) {
			return false;
		}
	} else {
		// Setup bones for current frame
		matrix3x4 aBones[MAXSTUDIOBONES];
		if (!tTarget.m_pEntity->SetupBones(aBones, MAXSTUDIOBONES, BONE_USED_BY_ANYTHING, tTarget.m_pEntity->m_flSimulationTime())) {
			return false;
		}

		F::Backtrack.m_tRecord = TickRecord{
			tTarget.m_pEntity->m_flSimulationTime(),
			tTarget.m_pEntity->m_vecOrigin(),
			tTarget.m_pEntity->m_vecMins(),
			tTarget.m_pEntity->m_vecMaxs(),
			*reinterpret_cast<BoneMatrix*>(&aBones),
			false,
			Vec3{},
			false
		};
		vRecords = { &F::Backtrack.m_tRecord };
	}

	CGameTrace trace{};
	CTraceFilterHitscan filter{};
	filter.pSkip = pLocal;
	
	const bool isKnife = pWeapon->GetWeaponID() == TF_WEAPON_KNIFE;
	const bool autoBackstab = Vars::Aimbot::Melee::AutoBackstab.Value;
	const int aimType = Vars::Aimbot::General::AimType.Value;
	constexpr float compressionOffset = 0.125f;
	
	for (const auto pRecord : vRecords) {
		// Store original values for restoration
		const Vec3 vRestoreOrigin = tTarget.m_pEntity->GetAbsOrigin();
		const Vec3 vRestoreMins = tTarget.m_pEntity->m_vecMins();
		const Vec3 vRestoreMaxs = tTarget.m_pEntity->m_vecMaxs();

		// Set entity to record position
		tTarget.m_pEntity->SetAbsOrigin(pRecord->m_vOrigin);
		tTarget.m_pEntity->m_vecMins() = pRecord->m_vMins + compressionOffset; // account for origin compression
		tTarget.m_pEntity->m_vecMaxs() = pRecord->m_vMaxs - compressionOffset;

		// Calculate target position
		const Vec3 vDiff{ 0, 0, std::clamp(vEyePos.z - pRecord->m_vOrigin.z,
										   tTarget.m_pEntity->m_vecMins().z,
										   tTarget.m_pEntity->m_vecMaxs().z) };
		tTarget.m_vPos = pRecord->m_vOrigin + vDiff;
		Aim(G::CurrentUserCmd->viewangles, Math::CalcAngle(vEyePos, tTarget.m_vPos), tTarget.m_vAngleTo);

		Vec3 vForward;
		Math::AngleVectors(tTarget.m_vAngleTo, &vForward);
		const Vec3 vTraceEnd = vEyePos + (vForward * flRange);

		// Primary trace
		SDK::TraceHull(vEyePos, vTraceEnd, {}, {}, MASK_SOLID, &filter, &trace);
		bool bReturn = trace.m_pEnt && trace.m_pEnt == tTarget.m_pEntity;
		
		// Secondary trace with swing bounds if primary failed
		if (!bReturn) {
			SDK::TraceHull(vEyePos, vTraceEnd, vSwingMins, vSwingMaxs, MASK_SOLID, &filter, &trace);
			bReturn = trace.m_pEnt && trace.m_pEnt == tTarget.m_pEntity;
		}

		// Check backstab conditions
		if (bReturn && autoBackstab && isKnife) {
			bReturn = (tTarget.m_iTargetType == TargetEnum::Player)
				? CanBackstab(tTarget.m_pEntity, pLocal, tTarget.m_vAngleTo)
				: false;
		}

		// Restore entity state
		tTarget.m_pEntity->SetAbsOrigin(vRestoreOrigin);
		tTarget.m_pEntity->m_vecMins() = vRestoreMins;
		tTarget.m_pEntity->m_vecMaxs() = vRestoreMaxs;

		if (bReturn) {
			tTarget.m_pRecord = pRecord;
			tTarget.m_bBacktrack = tTarget.m_iTargetType == TargetEnum::Player;
			return true;
		}
		
		// Handle smooth/assistive aim types
		if (aimType == Vars::Aimbot::General::AimTypeEnum::Smooth ||
			aimType == Vars::Aimbot::General::AimTypeEnum::Assistive) {
			const Vec3 vAngle = Math::CalcAngle(vEyePos, tTarget.m_vPos);
			Vec3 vForwardAlt;
			Math::AngleVectors(vAngle, &vForwardAlt);
			const Vec3 vTraceEndAlt = vEyePos + (vForwardAlt * flRange);

			SDK::Trace(vEyePos, vTraceEndAlt, MASK_SHOT | CONTENTS_GRATE, &filter, &trace);
			if (trace.m_pEnt && trace.m_pEnt == tTarget.m_pEntity) {
				return 2;
			}
		}
	}

	return false;
}



bool CAimbotMelee::Aim(Vec3 vCurAngle, Vec3 vToAngle, Vec3& vOut, int iMethod)
{
	if (const Vec3* pDoubletapAngle = F::Ticks.GetShootAngle()) {
		vOut = *pDoubletapAngle;
		return true;
	}

	Math::ClampAngles(vToAngle);

	switch (iMethod) {
	case Vars::Aimbot::General::AimTypeEnum::Plain:
	case Vars::Aimbot::General::AimTypeEnum::Silent:
	case Vars::Aimbot::General::AimTypeEnum::Locking:
		vOut = vToAngle;
		return false;
		
	case Vars::Aimbot::General::AimTypeEnum::Smooth: {
		const float assistStrength = Vars::Aimbot::General::AssistStrength.Value / 100.f;
		vOut = vCurAngle.LerpAngle(vToAngle, assistStrength);
		return true;
	}
	
	case Vars::Aimbot::General::AimTypeEnum::Assistive: {
		const Vec3 vMouseDelta = G::CurrentUserCmd->viewangles.DeltaAngle(G::LastUserCmd->viewangles);
		Vec3 vTargetDelta = vToAngle.DeltaAngle(G::LastUserCmd->viewangles);
		const float flMouseDelta = vMouseDelta.Length2D();
		const float flTargetDelta = vTargetDelta.Length2D();
		const Vec3 vNormalizedTargetDelta = vTargetDelta.Normalized() * std::min(flMouseDelta, flTargetDelta);
		const float assistStrength = Vars::Aimbot::General::AssistStrength.Value / 100.f;
		
		vOut = vCurAngle - vMouseDelta + vMouseDelta.LerpAngle(vNormalizedTargetDelta, assistStrength);
		return true;
	}
	
	default:
		return false;
	}
}

// assume angle calculated outside with other overload
void CAimbotMelee::Aim(CUserCmd* pCmd, Vec3& vAngle)
{
	switch (Vars::Aimbot::General::AimType.Value)
	{
	case Vars::Aimbot::General::AimTypeEnum::Plain:
	case Vars::Aimbot::General::AimTypeEnum::Smooth:
	case Vars::Aimbot::General::AimTypeEnum::Assistive:
		pCmd->viewangles = vAngle;
		I::EngineClient->SetViewAngles(vAngle);
		break;
	case Vars::Aimbot::General::AimTypeEnum::Silent:
	{
		bool bDoubleTap = F::Ticks.m_bDoubletap || F::Ticks.GetTicks(H::Entities.GetWeapon()) || F::Ticks.m_bSpeedhack;
		if (G::Attacking == 1 || bDoubleTap)
		{
			SDK::FixMovement(pCmd, vAngle);
			pCmd->viewangles = vAngle;
			G::PSilentAngles = true;
		}
		break;
	}
	case Vars::Aimbot::General::AimTypeEnum::Locking:
	{
		SDK::FixMovement(pCmd, vAngle);
		pCmd->viewangles = vAngle;
	}
	}
}

void CAimbotMelee::Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	static int iStaticAimType = Vars::Aimbot::General::AimType.Value;
	const int iLastAimType = iStaticAimType;
	const int iRealAimType = Vars::Aimbot::General::AimType.Value;

	// Handle aim type persistence during weapon swing
	if (pWeapon->m_flSmackTime() > 0.f && !iRealAimType && iLastAimType) {
		Vars::Aimbot::General::AimType.Value = iLastAimType;
	}
	iStaticAimType = Vars::Aimbot::General::AimType.Value;

	// Early returns for invalid conditions
	if (F::AimbotGlobal.ShouldHoldAttack(pWeapon)) {
		pCmd->buttons |= IN_ATTACK;
	}
	
	if (!Vars::Aimbot::General::AimType.Value ||
		(!F::AimbotGlobal.ShouldAim() && pWeapon->m_flSmackTime() < 0.f)) {
		return;
	}

	// Handle sapper logic
	if (RunSapper(pLocal, pWeapon, pCmd)) {
		return;
	}

	// Get and validate targets
	auto vTargets = SortTargets(pLocal, pWeapon);
	if (vTargets.empty()) {
		return;
	}

	// Calculate timing parameters
	iDoubletapTicks = F::Ticks.GetTicks(pWeapon);
	const int swingTime = GetSwingTime(pWeapon);
	const bool bShouldSwing = iDoubletapTicks <= (swingTime ? 14 : 0) ||
							  (Vars::Doubletap::AntiWarp.Value && pLocal->m_hGroundEntity());

	// Setup simulation
	Vec3 vEyePos = pLocal->GetShootPos();
	std::unordered_map<int, std::deque<TickRecord>> mRecordMap;
	std::unordered_map<int, std::vector<Vec3>> mPaths;
	mRecordMap.reserve(vTargets.size());
	mPaths.reserve(vTargets.size());
	
	SimulatePlayers(pLocal, pWeapon, vTargets, vEyePos, mRecordMap, mPaths);

	//if (!G::AimTarget.m_iEntIndex)
	//	G::AimTarget = { vTargets.front().m_pEntity->entindex(), I::GlobalVars->tickcount, 0 };

	for (auto& tTarget : vTargets)
	{
		auto iResult = CanHit(tTarget, pLocal, pWeapon, vEyePos, mRecordMap[tTarget.m_pEntity->entindex()]);
		if (!iResult) continue;
		if (iResult == 2)
		{
			G::AimTarget = { tTarget.m_pEntity->entindex(), I::GlobalVars->tickcount, 0 };
			Aim(pCmd, tTarget.m_vAngleTo);
			break;
		}

		G::AimTarget = { tTarget.m_pEntity->entindex(), I::GlobalVars->tickcount };
		G::AimPoint = { tTarget.m_vPos, I::GlobalVars->tickcount };

		if (Vars::Aimbot::General::AutoShoot.Value && pWeapon->m_flSmackTime() < 0.f)
		{
			if (bShouldSwing)
				pCmd->buttons |= IN_ATTACK;
			if (iDoubletapTicks)
				F::Ticks.m_bDoubletap = true;
		}

		G::Attacking = SDK::IsAttacking(pLocal, pWeapon, pCmd, true);

		if (G::Attacking == 1)
		{
			if (tTarget.m_bBacktrack)
				pCmd->tick_count = TIME_TO_TICKS(tTarget.m_pRecord->m_flSimTime) + TIME_TO_TICKS(F::Backtrack.GetFakeInterp());
			// bug: fast old records seem to be progressively more unreliable ?
		}
		else
		{
			const Vec3 vLocalEyePos = pLocal->GetShootPos();
			Aim(G::CurrentUserCmd->viewangles, Math::CalcAngle(vLocalEyePos, tTarget.m_vPos), tTarget.m_vAngleTo);
		}

		bool bPath = Vars::Visuals::Simulation::SwingLines.Value && Vars::Visuals::Simulation::PlayerPath.Value && Vars::Aimbot::General::AutoShoot.Value && !Vars::Debug::Info.Value;
		bool bLine = Vars::Visuals::Line::Enabled.Value;
		bool bBoxes = Vars::Visuals::Hitbox::BonesEnabled.Value & Vars::Visuals::Hitbox::BonesEnabledEnum::OnShot;
		if (pCmd->buttons & IN_ATTACK && pWeapon->m_flSmackTime() < 0.f && bPath)
		{
			G::LineStorage.clear();
			G::BoxStorage.clear();
			G::PathStorage.clear();

			if (bPath)
			{
				if (Vars::Colors::PlayerPath.Value.a)
				{
					G::PathStorage.emplace_back(mPaths[pLocal->entindex()], I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
					G::PathStorage.emplace_back(mPaths[tTarget.m_pEntity->entindex()], I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
				}
				if (Vars::Colors::PlayerPathClipped.Value.a)
				{
					G::PathStorage.emplace_back(mPaths[pLocal->entindex()], I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
					G::PathStorage.emplace_back(mPaths[tTarget.m_pEntity->entindex()], I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
				}
			}
		}
		if (G::Attacking == 1 && (bLine || bBoxes))
		{
			G::LineStorage.clear();
			G::BoxStorage.clear();

			if (bLine)
			{
				Vec3 vEyePos = pLocal->GetShootPos();
				float flDist = vEyePos.DistTo(tTarget.m_vPos);
				Vec3 vForward; Math::AngleVectors(tTarget.m_vAngleTo + pLocal->m_vecPunchAngle(), &vForward);

				if (Vars::Colors::Line.Value.a)
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(vEyePos, vEyePos + vForward * flDist), I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::Line.Value);
				if (Vars::Colors::LineClipped.Value.a)
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(vEyePos, vEyePos + vForward * flDist), I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::LineClipped.Value, true);
			}
			if (bBoxes)
			{
				auto vBoxes = F::Visuals.GetHitboxes(tTarget.m_pRecord->m_BoneMatrix.m_aBones, tTarget.m_pEntity->As<CBaseAnimating>());
				G::BoxStorage.insert(G::BoxStorage.end(), vBoxes.begin(), vBoxes.end());

				//if (Vars::Colors::BoneHitboxEdge.Value.a || Vars::Colors::BoneHitboxFace.Value.a)
				//	G::BoxStorage.emplace_back(tTarget.m_pRecord->m_vOrigin, tTarget.m_pRecord->m_vMins, tTarget.m_pRecord->m_vMaxs, Vec3(), I::GlobalVars->curtime + Vars::Visuals::Hitbox::DrawDuration.Value, Vars::Colors::BoneHitboxEdge.Value, Vars::Colors::BoneHitboxFace.Value);
				//if (Vars::Colors::BoneHitboxEdgeClipped.Value.a || Vars::Colors::BoneHitboxFaceClipped.Value.a)
				//	G::BoxStorage.emplace_back(tTarget.m_pRecord->m_vOrigin, tTarget.m_pRecord->m_vMins, tTarget.m_pRecord->m_vMaxs, Vec3(), I::GlobalVars->curtime + Vars::Visuals::Hitbox::DrawDuration.Value, Vars::Colors::BoneHitboxEdgeClipped.Value, Vars::Colors::BoneHitboxFaceClipped.Value, true);
			}
		}

		Aim(pCmd, tTarget.m_vAngleTo);
		break;
	}
}

static inline int GetAttachment(CBaseObject* pBuilding, int i)
{
	int iAttachment = pBuilding->GetBuildPointAttachmentIndex(i);
	if (pBuilding->IsSentrygun() && pBuilding->m_iUpgradeLevel() > 1) // idk why i need this
		iAttachment = 3;
	return iAttachment;
}
bool CAimbotMelee::FindNearestBuildPoint(CBaseObject* pBuilding, CTFPlayer* pLocal, Vec3& vPoint)
{
	bool bFoundPoint = false;

	static auto tf_obj_max_attach_dist = U::ConVars.FindVar("tf_obj_max_attach_dist");
	float flNearestPoint = tf_obj_max_attach_dist->GetFloat();
	for (int i = 0; i < pBuilding->GetNumBuildPoints(); i++)
	{
		int v = GetAttachment(pBuilding, i);

		Vec3 vOrigin;
		if (pBuilding->GetAttachment(v, vOrigin)) // issues using pBuilding->GetBuildPoint i on sentries above level 1 for some reason
		{
			if (!SDK::VisPos(pLocal, pBuilding, pLocal->GetShootPos(), vOrigin))
				continue;

			float flDist = (vOrigin - pLocal->GetAbsOrigin()).Length();
			if (flDist < flNearestPoint)
			{
				flNearestPoint = flDist;
				vPoint = vOrigin;
				bFoundPoint = true;
			}
		}
	}

	return bFoundPoint;
}

bool CAimbotMelee::RunSapper(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	if (pWeapon->GetWeaponID() != TF_WEAPON_BUILDER) {
		return false;
	}

	const Vec3 vLocalPos = pLocal->GetShootPos();
	const Vec3 vLocalAngles = I::EngineClient->GetViewAngles();
	const float maxFOV = Vars::Aimbot::General::AimFOV.Value;

	std::vector<Target_t> vTargets;
	vTargets.reserve(16); // Reserve space for enemy buildings
	
	const auto& enemyBuildings = H::Entities.GetGroup(EGroupType::BUILDINGS_ENEMIES);
	for (const auto pEntity : enemyBuildings) {
		const auto pBuilding = pEntity->As<CBaseObject>();
		if (pBuilding->m_bHasSapper() || !pBuilding->IsInValidTeam()) {
			continue;
		}

		Vec3 vPoint;
		if (!FindNearestBuildPoint(pBuilding, pLocal, vPoint)) {
			continue;
		}

		const Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPoint);
		const float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
		if (flFOVTo > maxFOV) {
			continue;
		}

		const float flDistTo = vLocalPos.DistTo(vPoint);
		vTargets.emplace_back(pBuilding, TargetEnum::Unknown, vPoint, vAngleTo, flFOVTo, flDistTo);
	}
	
	F::AimbotGlobal.SortTargets(vTargets, Vars::Aimbot::General::TargetSelectionEnum::Distance);
	if (vTargets.empty()) {
		return true;
	}

	//if (!G::AimTarget.m_iEntIndex)
	//	G::AimTarget = { vTargets.front().m_pEntity->entindex(), I::GlobalVars->tickcount, 0 };

	// Process first valid target
	for (auto& tTarget : vTargets) {
		static int iLastRun = 0;
		const int currentTick = I::GlobalVars->tickcount;
		const bool autoShoot = Vars::Aimbot::General::AutoShoot.Value;
		const int aimType = Vars::Aimbot::General::AimType.Value;

		bool bShouldAim = true;
		if (autoShoot) {
			pCmd->buttons |= IN_ATTACK;
		} else {
			bShouldAim = pCmd->buttons & IN_ATTACK;
		}
		
		if (aimType == Vars::Aimbot::General::AimTypeEnum::Silent) {
			bShouldAim = bShouldAim && (iLastRun != currentTick - 1 || (G::PSilentAngles && !F::Ticks.CanChoke()));
		}
		
		if (bShouldAim) {
			G::AimTarget = { tTarget.m_pEntity->entindex(), currentTick };
			G::AimPoint = { tTarget.m_vPos, currentTick };
			G::Attacking = true;

			Aim(pCmd->viewangles, Math::CalcAngle(vLocalPos, tTarget.m_vPos), tTarget.m_vAngleTo);
			tTarget.m_vAngleTo.x = pCmd->viewangles.x; // we don't need to care about pitch
			Aim(pCmd, tTarget.m_vAngleTo);

			iLastRun = currentTick;
		}

		break; // Only process first target
	}

	return true;
}