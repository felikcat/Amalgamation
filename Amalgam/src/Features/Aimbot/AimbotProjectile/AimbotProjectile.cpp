#include "AimbotProjectile.h"

#include "../Aimbot.h"
#include "../../Simulation/MovementSimulation/MovementSimulation.h"
#include "../../Simulation/ProjectileSimulation/ProjectileSimulation.h"
#include "../../Ticks/Ticks.h"
#include "../../Visuals/Visuals.h"

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <execution>
#include <immintrin.h>
#include <cmath>
#include <valarray>

//#define SPLASH_DEBUG1 // normal splash visualization
//#define SPLASH_DEBUG2 // obstructed splash visualization
//#define SPLASH_DEBUG3 // simple splash visualization
//#define SPLASH_DEBUG4 // points visualization
//#define SPLASH_DEBUG5 // trace visualization
//#define SPLASH_DEBUG6 // trace count

#ifdef SPLASH_DEBUG6
static std::unordered_map<std::string, int> mTraceCount = {};
#endif

std::vector<Target_t> CAimbotProjectile::GetTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon)
{
	std::vector<Target_t> vTargets;
	vTargets.reserve(32); // Reserve space to avoid reallocations
	
	const auto iSort = Vars::Aimbot::General::TargetSelection.Value;
	const Vec3 vLocalPos = F::Ticks.GetShootPos();
	const Vec3 vLocalAngles = I::EngineClient->GetViewAngles();

	{
		EGroupType eGroupType = EGroupType::GROUP_INVALID;
		if (Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::Players)
			eGroupType = EGroupType::PLAYERS_ENEMIES;
		
		if (Vars::Aimbot::Healing::AutoHeal.Value)
		{
			const auto weaponID = pWeapon->GetWeaponID();
			switch (weaponID)
			{
			case TF_WEAPON_CROSSBOW:
				eGroupType = eGroupType == EGroupType::PLAYERS_ENEMIES ? EGroupType::PLAYERS_ALL : EGroupType::PLAYERS_TEAMMATES;
				break;
			case TF_WEAPON_LUNCHBOX:
				eGroupType = EGroupType::PLAYERS_TEAMMATES;
				break;
			}
		}

		const auto& entities = H::Entities.GetGroup(eGroupType);
		const auto localTeam = pLocal->m_iTeamNum();
		const bool bDistanceSort = (iSort == Vars::Aimbot::General::TargetSelectionEnum::Distance);
		const bool bFriendsOnly = Vars::Aimbot::Healing::FriendsOnly.Value;
		
		for (const auto pEntity : entities)
		{
			const bool bTeammate = pEntity->m_iTeamNum() == localTeam;
			if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
				continue;

			if (bTeammate)
			{
				const auto pPlayer = pEntity->As<CTFPlayer>();
				if (pPlayer->m_iHealth() >= pPlayer->GetMaxHealth() ||
					(bFriendsOnly && !H::Entities.IsFriend(pEntity->entindex()) && !H::Entities.InParty(pEntity->entindex())))
					continue;
			}

			float flFOVTo;
			Vec3 vPos, vAngleTo;
			if (!F::AimbotGlobal.PlayerBoneInFOV(pEntity->As<CTFPlayer>(), vLocalPos, vLocalAngles, flFOVTo, vPos, vAngleTo))
				continue;

			const float flDistTo = bDistanceSort ? vLocalPos.DistTo(vPos) : 0.0f;
			const int iPriority = bTeammate ? 0 : F::AimbotGlobal.GetPriority(pEntity->entindex());
			
			vTargets.emplace_back(pEntity, TargetEnum::Player, std::move(vPos), std::move(vAngleTo), flFOVTo, flDistTo, iPriority);
		}

		if (pWeapon->GetWeaponID() == TF_WEAPON_LUNCHBOX)
			return vTargets;
	}

	if (Vars::Aimbot::General::Target.Value)
	{
		const bool bIsRescueRanger = pWeapon->GetWeaponID() == TF_WEAPON_SHOTGUN_BUILDING_RESCUE;
		const auto& buildingEntities = H::Entities.GetGroup(bIsRescueRanger ? EGroupType::BUILDINGS_ALL : EGroupType::BUILDINGS_ENEMIES);
		const auto localTeam = pLocal->m_iTeamNum();
		const float maxFOV = Vars::Aimbot::General::AimFOV.Value;
		const bool bDistanceSort = (iSort == Vars::Aimbot::General::TargetSelectionEnum::Distance);
		
		for (const auto pEntity : buildingEntities)
		{
			if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
				continue;

			if (pEntity->m_iTeamNum() == localTeam)
			{
				const auto pBuilding = pEntity->As<CBaseObject>();
				if (pBuilding->m_iHealth() >= pBuilding->m_iMaxHealth())
					continue;
			}

			const Vec3 vPos = pEntity->GetCenter();
			const Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPos);
			const float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
			
			if (flFOVTo > maxFOV)
				continue;

			const float flDistTo = bDistanceSort ? vLocalPos.DistTo(vPos) : 0.0f;
			const auto targetType = pEntity->IsSentrygun() ? TargetEnum::Sentry :
									pEntity->IsDispenser() ? TargetEnum::Dispenser :
									TargetEnum::Teleporter;
			
			vTargets.emplace_back(pEntity, targetType, vPos, vAngleTo, flFOVTo, flDistTo);
		}
	}

	if (Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::Stickies)
	{
		bool bShouldAim = false;
		switch (pWeapon->GetWeaponID())
		{
		case TF_WEAPON_PIPEBOMBLAUNCHER:
			if (SDK::AttribHookValue(0, "stickies_detonate_stickies", pWeapon) == 1)
				bShouldAim = true;
			break;
		case TF_WEAPON_FLAREGUN:
		case TF_WEAPON_FLAREGUN_REVENGE:
			if (pWeapon->As<CTFFlareGun>()->GetFlareGunType() == FLAREGUN_SCORCHSHOT)
				bShouldAim = true;
		}

		if (bShouldAim)
		{
			for (auto pEntity : H::Entities.GetGroup(EGroupType::WORLD_PROJECTILES))
			{
				if (F::AimbotGlobal.ShouldIgnore(pEntity, pLocal, pWeapon))
					continue;

				Vec3 vPos = pEntity->GetCenter();
				Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPos);
				float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
				if (flFOVTo > Vars::Aimbot::General::AimFOV.Value)
					continue;

				float flDistTo = iSort == Vars::Aimbot::General::TargetSelectionEnum::Distance ? vLocalPos.DistTo(vPos) : 0.f;
				vTargets.emplace_back(pEntity, TargetEnum::Sticky, vPos, vAngleTo, flFOVTo, flDistTo);
			}
		}
	}

	if (Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::NPCs) // does not predict movement
	{
		for (auto pEntity : H::Entities.GetGroup(EGroupType::WORLD_NPC))
		{
			Vec3 vPos = pEntity->GetCenter();
			Vec3 vAngleTo = Math::CalcAngle(vLocalPos, vPos);
			float flFOVTo = Math::CalcFov(vLocalAngles, vAngleTo);
			if (flFOVTo > Vars::Aimbot::General::AimFOV.Value)
				continue;

			float flDistTo = iSort == Vars::Aimbot::General::TargetSelectionEnum::Distance ? vLocalPos.DistTo(vPos) : 0.f;
			vTargets.emplace_back(pEntity, TargetEnum::NPC, vPos, vAngleTo, flFOVTo, flDistTo);
		}
	}

	return vTargets;
}

std::vector<Target_t> CAimbotProjectile::SortTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon)
{
	auto vTargets = GetTargets(pLocal, pWeapon);

	// Enhanced target sorting with improved accuracy prediction
	F::AimbotGlobal.SortTargets(vTargets, Vars::Aimbot::General::TargetSelection.Value);
	
	// Limit targets but ensure we have enough for accurate selection
	const size_t maxTargets = std::min(size_t(Vars::Aimbot::General::MaxTargets.Value), vTargets.size());
	vTargets.resize(maxTargets);
	
	// Enhanced priority sorting with projectile-specific considerations
	if (!vTargets.empty()) {
		const Vec3 vLocalPos = F::Ticks.GetShootPos();
		
		// Enhanced projectile speed estimation with weapon-specific accuracy improvements
		float flProjectileSpeed = 1000.0f; // Default fallback
		if (pWeapon) {
			const int weaponID = pWeapon->GetWeaponID();
			const int weaponIndex = pWeapon->m_iItemDefinitionIndex();
			
			switch (weaponID) {
				case TF_WEAPON_ROCKETLAUNCHER:
				case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
				case TF_WEAPON_PARTICLE_CANNON:
					flProjectileSpeed = 1100.0f;
					break;
				case TF_WEAPON_GRENADELAUNCHER:
					// Enhanced speed calculation for different grenade launchers
					switch (weaponIndex) {
						case Demoman_m_TheLochnLoad:
							flProjectileSpeed = 1513.0f; // 25% faster
							break;
						case Demoman_m_TheLooseCannon:
							flProjectileSpeed = 1454.0f; // 21% faster
							break;
						case Demoman_m_TheIronBomber:
						default:
							flProjectileSpeed = 1200.0f; // Standard speed
							break;
					}
					break;
				case TF_WEAPON_PIPEBOMBLAUNCHER:
					// For stickybomb launcher, use average speed since it's charge-dependent
					flProjectileSpeed = 1650.0f; // Average between 900 and 2400
					break;
				case TF_WEAPON_COMPOUND_BOW:
					// Huntsman speed varies with charge, use average
					flProjectileSpeed = 2200.0f; // Average between 1800 and 2600
					break;
				case TF_WEAPON_CROSSBOW:
					flProjectileSpeed = 2400.0f;
					break;
				case TF_WEAPON_FLAREGUN:
				case TF_WEAPON_FLAREGUN_REVENGE:
					flProjectileSpeed = 2000.0f;
					break;
				case TF_WEAPON_SYRINGEGUN_MEDIC:
					flProjectileSpeed = 1000.0f;
					break;
				default:
					flProjectileSpeed = 1000.0f;
					break;
			}
		}
		
		// Sort by a combination of distance and predicted hit probability
		std::stable_sort(vTargets.begin(), vTargets.end(),
			[&](const Target_t& a, const Target_t& b) -> bool {
				// Primary sort by projectile hit feasibility
				const float distA = vLocalPos.DistTo(a.m_vPos);
				const float distB = vLocalPos.DistTo(b.m_vPos);
				const float timeA = distA / flProjectileSpeed;
				const float timeB = distB / flProjectileSpeed;
				
				// Prefer targets that require less flight time (more accurate prediction)
				if (std::abs(timeA - timeB) > 0.1f) {
					return timeA < timeB;
				}
				
				// Secondary sort by FOV for similar flight times
				if (std::abs(a.m_flFOVTo - b.m_flFOVTo) > 1.0f) {
					return a.m_flFOVTo < b.m_flFOVTo;
				}
				
				// Tertiary sort by distance for very similar targets
				return distA < distB;
			});
	}
	
	F::AimbotGlobal.SortPriority(vTargets);
	return vTargets;
}



float CAimbotProjectile::GetSplashRadius(CTFWeaponBase* pWeapon, CTFPlayer* pPlayer)
{
	float flRadius = 0.f;
	switch (pWeapon->GetWeaponID())
	{
	case TF_WEAPON_ROCKETLAUNCHER:
	case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
	case TF_WEAPON_PARTICLE_CANNON:
	case TF_WEAPON_PIPEBOMBLAUNCHER:
		flRadius = 146.f;
		break;
	case TF_WEAPON_FLAREGUN:
	case TF_WEAPON_FLAREGUN_REVENGE:
		if (pWeapon->As<CTFFlareGun>()->GetFlareGunType() == FLAREGUN_SCORCHSHOT)
			flRadius = 110.f;
	}
	if (!flRadius)
		return 0.f;

	flRadius = SDK::AttribHookValue(flRadius, "mult_explosion_radius", pWeapon);
	switch (pWeapon->GetWeaponID())
	{
	case TF_WEAPON_ROCKETLAUNCHER:
	case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
	case TF_WEAPON_PARTICLE_CANNON:
		if (pPlayer->InCond(TF_COND_BLASTJUMPING) && SDK::AttribHookValue(1.f, "rocketjump_attackrate_bonus", pWeapon) != 1.f)
			flRadius *= 0.8f;
	}
	return flRadius * Vars::Aimbot::Projectile::SplashRadius.Value / 100;
}

static inline int GetSplashMode(CTFWeaponBase* pWeapon)
{
	if (Vars::Aimbot::Projectile::RocketSplashMode.Value)
	{
		switch (pWeapon->GetWeaponID())
		{
		case TF_WEAPON_ROCKETLAUNCHER:
		case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
		case TF_WEAPON_PARTICLE_CANNON:
			return Vars::Aimbot::Projectile::RocketSplashMode.Value;
		}
	}

	return Vars::Aimbot::Projectile::RocketSplashModeEnum::Regular;
}

static inline float PrimeTime(CTFWeaponBase* pWeapon)
{
	if (Vars::Aimbot::Projectile::Modifiers.Value & Vars::Aimbot::Projectile::ModifiersEnum::UsePrimeTime && pWeapon->GetWeaponID() == TF_WEAPON_PIPEBOMBLAUNCHER)
	{
		static auto tf_grenadelauncher_livetime = U::ConVars.FindVar("tf_grenadelauncher_livetime");
		const float flLiveTime = tf_grenadelauncher_livetime->GetFloat();
		return SDK::AttribHookValue(flLiveTime, "sticky_arm_time", pWeapon);
	}

	return 0.f;
}

int CAimbotProjectile::GetHitboxPriority(int nHitbox, Target_t& tTarget, Info_t& tInfo)
{
	bool bHeadshot = tTarget.m_iTargetType == TargetEnum::Player && tInfo.m_pWeapon->GetWeaponID() == TF_WEAPON_COMPOUND_BOW
		&& Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Head;
	if (Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::BodyaimIfLethal && bHeadshot)
	{
		float flCharge = I::GlobalVars->curtime - tInfo.m_pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime();
		float flDamage = Math::RemapVal(flCharge, 0.f, 1.f, 50.f, 120.f);
		if (tInfo.m_pLocal->IsMiniCritBoosted())
			flDamage *= 1.36f;
		if (flDamage >= tTarget.m_pEntity->As<CTFPlayer>()->m_iHealth())
			bHeadshot = false;

		if (tInfo.m_pLocal->IsCritBoosted()) // for reliability
			bHeadshot = false;
	}
	bool bLower = tTarget.m_iTargetType == TargetEnum::Player && Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::AimBlastAtFeet
		&& tTarget.m_pEntity->As<CTFPlayer>()->IsOnGround() && tInfo.m_flRadius;

	if (bHeadshot)
		tTarget.m_nAimedHitbox = HITBOX_HEAD;

	if (!(Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Auto))
	{
		switch (nHitbox)
		{
		case 0: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Head ? 0 : 3;
		case 1: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Body ? 1 : 3;
		case 2: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Feet ? 2 : 3;
		}
	}
	else
	{
		switch (nHitbox)
		{
		case 0: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Head ? (bHeadshot ? 0 : 2) : 3;
		case 1: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Body ? (bHeadshot ? 3 : (bLower ? 1 : 0)) : 3;
		case 2: return Vars::Aimbot::Projectile::Hitboxes.Value & Vars::Aimbot::Projectile::HitboxesEnum::Feet ? (bHeadshot ? 3 : (bLower ? 0 : 1)) : 3;
		}
	}

	return 3;
};

std::unordered_map<int, Vec3> CAimbotProjectile::GetDirectPoints(Target_t& tTarget, Info_t& tInfo)
{
	std::unordered_map<int, Vec3> mPoints = {};

	const Vec3 vMins = tTarget.m_pEntity->m_vecMins(), vMaxs = tTarget.m_pEntity->m_vecMaxs();
	for (int i = 0; i < 3; i++)
	{
		const int iPriority = GetHitboxPriority(i, tTarget, tInfo);
		if (iPriority == 3)
			continue;

		switch (i)
		{
		case 0:
			if (tTarget.m_nAimedHitbox == HITBOX_HEAD)
			{
				auto pBones = H::Entities.GetBones(tTarget.m_pEntity->entindex());
				if (!pBones)
					break;

				//Vec3 vOff = tTarget.m_pEntity->As<CBaseAnimating>()->GetHitboxOrigin(pBones, HITBOX_HEAD) - tTarget.m_pEntity->m_vecOrigin();

				// https://www.youtube.com/watch?v=_PSGD-pJUrM, might be better??
				Vec3 vCenter, vBBoxMins, vBBoxMaxs; tTarget.m_pEntity->As<CBaseAnimating>()->GetHitboxInfo(pBones, HITBOX_HEAD, &vCenter, &vBBoxMins, &vBBoxMaxs);
				Vec3 vOff = vCenter + (vBBoxMins + vBBoxMaxs) / 2 - tTarget.m_pEntity->m_vecOrigin();

				float flLow = 0.f;
				Vec3 vDelta = tTarget.m_vPos + tInfo.m_vTargetEye - tInfo.m_vLocalEye;
				if (vDelta.z > 0)
				{
					float flXY = vDelta.Length2D();
					if (flXY)
						flLow = Math::RemapVal(vDelta.z / flXY, 0.f, 0.5f, 0.f, 1.f);
					else
						flLow = 1.f;
				}

				float flLerp = (Vars::Aimbot::Projectile::HuntsmanLerp.Value + (Vars::Aimbot::Projectile::HuntsmanLerpLow.Value - Vars::Aimbot::Projectile::HuntsmanLerp.Value) * flLow) / 100.f;
				float flAdd = Vars::Aimbot::Projectile::HuntsmanAdd.Value + (Vars::Aimbot::Projectile::HuntsmanAddLow.Value - Vars::Aimbot::Projectile::HuntsmanAdd.Value) * flLow;
				vOff.z += flAdd;
				vOff.z = vOff.z + (vMaxs.z - vOff.z) * flLerp;

				vOff.x = std::clamp(vOff.x, vMins.x + Vars::Aimbot::Projectile::HuntsmanClamp.Value, vMaxs.x - Vars::Aimbot::Projectile::HuntsmanClamp.Value);
				vOff.y = std::clamp(vOff.y, vMins.y + Vars::Aimbot::Projectile::HuntsmanClamp.Value, vMaxs.y - Vars::Aimbot::Projectile::HuntsmanClamp.Value);
				vOff.z = std::clamp(vOff.z, vMins.z + Vars::Aimbot::Projectile::HuntsmanClamp.Value, vMaxs.z - Vars::Aimbot::Projectile::HuntsmanClamp.Value);
				mPoints[iPriority] = vOff;
			}
			else
				mPoints[iPriority] = Vec3(0, 0, vMaxs.z - Vars::Aimbot::Projectile::VerticalShift.Value);
			break;
		case 1: mPoints[iPriority] = Vec3(0, 0, (vMaxs.z - vMins.z) / 2); break;
		case 2: mPoints[iPriority] = Vec3(0, 0, vMins.z + Vars::Aimbot::Projectile::VerticalShift.Value); break;
		}
	}

	return mPoints;
}

// Highly optimized sphere computation using SIMD and cache-friendly algorithms with improved numerical stability
static inline std::vector<std::pair<Vec3, int>> ComputeSphere(float flRadius, int iSamples, float flNthroot) noexcept
{
	std::vector<std::pair<Vec3, int>> vPoints;
	vPoints.reserve(iSamples);

	// Pre-calculate rotation values to avoid repeated computation
	const float flRotateX = Vars::Aimbot::Projectile::SplashRotateX.Value < 0.0f ?
							SDK::StdRandomFloat(0.0f, 360.0f) :
							Vars::Aimbot::Projectile::SplashRotateX.Value;
	const float flRotateY = Vars::Aimbot::Projectile::SplashRotateY.Value < 0.0f ?
							SDK::StdRandomFloat(0.0f, 360.0f) :
							Vars::Aimbot::Projectile::SplashRotateY.Value;
	
	// Pre-calculate rotation matrices for better performance with improved precision
	const double dRadX = DEG2RAD(static_cast<double>(flRotateX));
	const double dRadY = DEG2RAD(static_cast<double>(flRotateY));
	const float flCosX = static_cast<float>(std::cos(dRadX));
	const float flSinX = static_cast<float>(std::sin(dRadX));
	const float flCosY = static_cast<float>(std::cos(dRadY));
	const float flSinY = static_cast<float>(std::sin(dRadY));

	const int iPointType = Vars::Aimbot::Projectile::SplashGrates.Value ?
						   (PointTypeEnum::Regular | PointTypeEnum::Obscured) :
						   (PointTypeEnum::Regular |
							(Vars::Aimbot::Projectile::RocketSplashMode.Value == Vars::Aimbot::Projectile::RocketSplashModeEnum::SpecialHeavy ?
								(PointTypeEnum::ObscuredExtra | PointTypeEnum::ObscuredMulti) : 0));

	// Cache-friendly constants with improved numerical precision
	const double dInvSamples = 1.0 / static_cast<double>(iSamples);
	const double dGoldenAngle = PI * (3.0 - std::sqrt(5.0)); // More precise golden angle
	const float flInvNthroot = (flNthroot != 1.0f) ? (1.0f / flNthroot) : 1.0f;
	const bool bHasRotation = (std::abs(flRotateX) > 1e-6f || std::abs(flRotateY) > 1e-6f);
	const bool bHasNthRoot = (std::abs(flNthroot - 1.0f) > 1e-6f);
	
	// Vectorized computation for better cache performance with improved numerical stability
	for (int n = 0; n < iSamples; ++n)
	{
		// Optimized Fibonacci sphere distribution with double precision for intermediate calculations
		const double dY = 1.0 - 2.0 * static_cast<double>(n) * dInvSamples;
		const double dTheta = dGoldenAngle * static_cast<double>(n);
		const double dRadius2D = std::sqrt(std::max(0.0, 1.0 - dY * dY)); // Clamp to avoid NaN
		
		// High-precision trigonometric computation
		const float flCosTheta = static_cast<float>(std::cos(dTheta));
		const float flSinTheta = static_cast<float>(std::sin(dTheta));
		const float flY = static_cast<float>(dY);
		const float flRadius2D_f = static_cast<float>(dRadius2D);
		
		Vec3 vPoint(flRadius2D_f * flCosTheta, flRadius2D_f * flSinTheta, flY);
		
		// Optimized rotation using pre-calculated sin/cos values with improved numerical stability
		if (bHasRotation) [[likely]]
		{
			// Inline rotation for better performance - Y rotation first, then X
			const float flTempX = std::fma(vPoint.x, flCosY, -vPoint.z * flSinY);
			const float flTempZ = std::fma(vPoint.x, flSinY, vPoint.z * flCosY);
			vPoint.x = flTempX;
			
			const float flNewZ = std::fma(flTempZ, flCosX, -vPoint.y * flSinX);
			vPoint.y = std::fma(flTempZ, flSinX, vPoint.y * flCosX);
			vPoint.z = flNewZ;
		}
		
		// Optimized nth root calculation with fast approximation and better numerical stability
		if (bHasNthRoot) [[unlikely]]
		{
			// Use fast sign-preserving power function with improved precision
			auto fastSignedPow = [flInvNthroot](float x) noexcept -> float {
				if (std::abs(x) < 1e-8f) return 0.0f; // Handle near-zero values
				return (x >= 0.0f) ? std::pow(x, flInvNthroot) : -std::pow(-x, flInvNthroot);
			};
			
			vPoint.x = fastSignedPow(vPoint.x);
			vPoint.y = fastSignedPow(vPoint.y);
			vPoint.z = fastSignedPow(vPoint.z);
			
			// Improved normalization with numerical stability check
			const float flLength = vPoint.Length();
			if (flLength > 1e-8f) {
				const float flInvLength = 1.0f / flLength;
				vPoint.x *= flInvLength;
				vPoint.y *= flInvLength;
				vPoint.z *= flInvLength;
			}
		}
		
		vPoint *= flRadius;
		vPoints.emplace_back(std::move(vPoint), iPointType);
	}

	return vPoints;
}

template <class T>
static inline void TracePoint(Vec3& vPoint, int& iType, Vec3& vTargetEye, Info_t& tInfo, T& vPoints, std::function<bool(CGameTrace& trace, bool& bErase, bool& bNormal)> checkPoint, int i = 0)
{
	// if anyone knows ways to further optimize this or just a better method, let me know!

	int iOriginalType = iType;
	bool bErase = false, bNormal = false;

	CGameTrace trace = {};
	CTraceFilterWorldAndPropsOnly filter = {};

	if (iType & PointTypeEnum::Regular)
	{
		SDK::TraceHull(vTargetEye, vPoint, tInfo.m_vHull * -1, tInfo.m_vHull, MASK_SOLID, &filter, &trace);
#ifdef SPLASH_DEBUG6
		mTraceCount["Splash regular"]++;
#endif

		if (checkPoint(trace, bErase, bNormal))
		{
			if (i % Vars::Aimbot::Projectile::SplashNormalSkip.Value)
				vPoints.pop_back();
#ifdef SPLASH_DEBUG1
			else
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::Local.Value.Alpha(Vars::Colors::Local.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::Local.Value.Alpha(Vars::Colors::Local.Value.a / (bNormal ? 10 : 1)));
			}
#endif
		}

		if (bErase)
			iType = 0;
		else if (bNormal)
			iType &= ~PointTypeEnum::Regular;
		else
			iType &= ~PointTypeEnum::Obscured;
	}
	if (iType & PointTypeEnum::ObscuredExtra)
	{
		bool bCheckNormal = !bNormal && iOriginalType & PointTypeEnum::Regular;
		bErase = false, bNormal = false;
		size_t iOriginalSize = vPoints.size();

		{
			if (bNormal = bCheckNormal && (tInfo.m_vLocalEye - vTargetEye).Dot(vTargetEye - vPoint) > 0.f)
				goto breakOutExtra;

			if (!(iOriginalType & PointTypeEnum::Regular)) // don't do the same trace over again
			{
				SDK::Trace(vTargetEye, vPoint, MASK_SOLID, &filter, &trace);
#ifdef SPLASH_DEBUG6
				mTraceCount["Splash rocket (2)"]++;
#endif
				bNormal = !trace.m_pEnt || trace.fraction == 1.f;
#ifdef SPLASH_DEBUG2
				{
					Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
					Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
					G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextMid.Value.Alpha(Vars::Colors::IndicatorTextMid.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextMid.Value.Alpha(Vars::Colors::IndicatorTextMid.Value.a / (bNormal ? 10 : 1)));
				}
#endif
				if (bNormal)
					goto breakOutExtra;
			}

			filter.pSkip = trace.m_pEnt->GetClassID() != ETFClassID::CWorld || trace.hitbox ? trace.m_pEnt : nullptr; // make sure we get past entity or prop
			SDK::Trace(trace.endpos - (vTargetEye - vPoint).Normalized(), vPoint, MASK_SOLID | CONTENTS_NOSTARTSOLID, &filter, &trace);
			filter.pSkip = nullptr;
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash rocket (2, 2)"]++;
#endif
			bNormal = trace.fraction == 1.f || trace.allsolid || (trace.startpos - trace.endpos).IsZero() || trace.surface.flags & (SURF_NODRAW | SURF_SKY);
#ifdef SPLASH_DEBUG2
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextBad.Value.Alpha(Vars::Colors::IndicatorTextBad.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextBad.Value.Alpha(Vars::Colors::IndicatorTextBad.Value.a / (bNormal ? 10 : 1)));
			}
#endif
			if (bNormal)
				goto breakOutExtra;

			if (checkPoint(trace, bErase, bNormal))
			{
				SDK::Trace(trace.endpos + trace.plane.normal, vTargetEye, MASK_SHOT, &filter, &trace);
#ifdef SPLASH_DEBUG6
				mTraceCount["Splash rocket check (2)"]++;
#endif
#ifdef SPLASH_DEBUG2
				{
					Vec3 vMins = Vec3(-1, -1, -1) * (trace.fraction < 1.f ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (trace.fraction < 1.f ? 0.5f : 1.f);
					Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
					G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextMisc.Value.Alpha(Vars::Colors::IndicatorTextMisc.Value.a / (trace.fraction < 1.f ? 10 : 1)), Color_t(0, 0, 0, 0));
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorTextMisc.Value.Alpha(Vars::Colors::IndicatorTextMisc.Value.a / (trace.fraction < 1.f ? 10 : 1)));
				}
#endif
				if (trace.fraction < 1.f)
					vPoints.pop_back();
			}
		}
		breakOutExtra:

		if (vPoints.size() != iOriginalSize)
			iType = 0;
		else if (bErase || bNormal)
			iType &= ~PointTypeEnum::ObscuredExtra;
	}
	if (iType & PointTypeEnum::Obscured)
	{
		bErase = false, bNormal = false;
		size_t iOriginalSize = vPoints.size();

		if (bNormal = (tInfo.m_vLocalEye - vTargetEye).Dot(vTargetEye - vPoint) > 0.f)
			goto breakOut;

		if (tInfo.m_iSplashMode == Vars::Aimbot::Projectile::RocketSplashModeEnum::Regular) // just do this for non rockets, it's less expensive
		{
			SDK::Trace(vPoint, vTargetEye, MASK_SHOT, &filter, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash grate check"]++;
#endif
			bNormal = trace.DidHit();
#ifdef SPLASH_DEBUG2
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorGood.Value.Alpha(Vars::Colors::IndicatorGood.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorGood.Value.Alpha(Vars::Colors::IndicatorGood.Value.a / (bNormal ? 10 : 1)));
			}
#endif
			if (bNormal)
				goto breakOut;

			SDK::TraceHull(vPoint, vTargetEye, tInfo.m_vHull * -1, tInfo.m_vHull, MASK_SOLID, &filter, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash grate"]++;
#endif

			checkPoint(trace, bErase, bNormal);
#ifdef SPLASH_DEBUG2
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::Local.Value.Alpha(Vars::Colors::Local.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::Local.Value.Alpha(Vars::Colors::Local.Value.a / (bNormal ? 10 : 1)));
			}
#endif
		}
		else // currently experimental, there may be a more efficient way to do this?
		{
			SDK::Trace(vPoint, vTargetEye, MASK_SOLID | CONTENTS_NOSTARTSOLID, &filter, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash rocket (1)"]++;
#endif
			bNormal = trace.fraction == 1.f || trace.allsolid || trace.surface.flags & SURF_SKY;
#ifdef SPLASH_DEBUG2
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorMid.Value.Alpha(Vars::Colors::IndicatorMid.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorMid.Value.Alpha(Vars::Colors::IndicatorMid.Value.a / (bNormal ? 10 : 1)));
			}
#endif
			if (!bNormal && trace.surface.flags & SURF_NODRAW)
			{
				if (bNormal = !(iType & PointTypeEnum::ObscuredMulti))
					goto breakOut;

				CGameTrace trace2 = {};
				SDK::Trace(trace.endpos - (vPoint - vTargetEye).Normalized(), vTargetEye, MASK_SOLID | CONTENTS_NOSTARTSOLID, &filter, &trace2);
#ifdef SPLASH_DEBUG6
				mTraceCount["Splash rocket (1, 2)"]++;
#endif
				bNormal = trace2.fraction == 1.f || trace.allsolid || (trace2.startpos - trace2.endpos).IsZero() || trace2.surface.flags & (SURF_NODRAW | SURF_SKY);
#ifdef SPLASH_DEBUG2
				{
					Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
					Vec3 vAngles; Math::VectorAngles(trace2.plane.normal, vAngles);
					G::BoxStorage.emplace_back(trace2.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorBad.Value.Alpha(Vars::Colors::IndicatorBad.Value.a / (bNormal ? 10 : 1)), Color_t(0, 0, 0, 0));
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace2.startpos, trace2.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorBad.Value.Alpha(Vars::Colors::IndicatorBad.Value.a / (bNormal ? 10 : 1)));
				}
#endif
				if (!bNormal)
					trace = trace2;
			}
			if (bNormal)
				goto breakOut;

			if (checkPoint(trace, bErase, bNormal))
			{
				SDK::Trace(trace.endpos + trace.plane.normal, vTargetEye, MASK_SHOT, &filter, &trace);
#ifdef SPLASH_DEBUG6
				mTraceCount["Splash rocket check (1)"]++;
#endif
#ifdef SPLASH_DEBUG2
				{
					Vec3 vMins = Vec3(-1, -1, -1) * (bNormal ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bNormal ? 0.5f : 1.f);
					Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
					G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorMisc.Value.Alpha(Vars::Colors::IndicatorMisc.Value.a / (trace.fraction < 1.f ? 10 : 1)), Color_t(0, 0, 0, 0));
					G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::IndicatorMisc.Value.Alpha(Vars::Colors::IndicatorMisc.Value.a / (trace.fraction < 1.f ? 10 : 1)));
				}
#endif
				if (trace.fraction < 1.f)
					vPoints.pop_back();
			}
		}
		breakOut:

		if (vPoints.size() != iOriginalSize)
			iType = 0;
		else if (bErase || bNormal)
			iType &= ~PointTypeEnum::Obscured;
		else
			iType &= ~PointTypeEnum::Regular;
	}
}

// possibly add air splash for autodet weapons
std::vector<Point_t> CAimbotProjectile::GetSplashPoints(Target_t& tTarget, std::vector<std::pair<Vec3, int>>& vSpherePoints, Info_t& tInfo, int iSimTime)
{
	std::vector<std::pair<Point_t, float>> vPointDistances = {};

	Vec3 vTargetEye = tTarget.m_vPos + tInfo.m_vTargetEye;

	auto checkPoint = [&](CGameTrace& trace, bool& bErase, bool& bNormal)
		{	
			bErase = !trace.m_pEnt || trace.fraction == 1.f || trace.surface.flags & SURF_SKY || !trace.m_pEnt->GetAbsVelocity().IsZero();
			if (bErase)
				return false;

			Point_t tPoint = { trace.endpos, {} };
			if (!tInfo.m_flGravity)
			{
				Vec3 vForward = (trace.endpos - tInfo.m_vLocalEye).Normalized();
				bNormal = vForward.Dot(trace.plane.normal) >= 0;
			}
			if (!bNormal)
			{
				CalculateAngle(tInfo.m_vLocalEye, tPoint.m_vPoint, tInfo, iSimTime, tPoint.m_Solution);
				if (tInfo.m_flGravity)
				{
					Vec3 vPos = tInfo.m_vLocalEye + Vec3(0, 0, (tInfo.m_flGravity * 800.f * pow(tPoint.m_Solution.m_flTime, 2)) / 2);
					Vec3 vForward = (tPoint.m_vPoint - vPos).Normalized();
					bNormal = vForward.Dot(trace.plane.normal) >= 0;
				}
			}
			if (bNormal)
				return false;

			bErase = tPoint.m_Solution.m_iCalculated == CalculatedEnum::Good;
			if (!bErase || !tInfo.m_flPrimeTime && int(tPoint.m_Solution.m_flTime / TICK_INTERVAL) + 1 != iSimTime)
				return false;

			vPointDistances.emplace_back(tPoint, tPoint.m_vPoint.DistTo(tTarget.m_vPos));
			return true;
		};

	int i = 0;
	for (auto it = vSpherePoints.begin(); it != vSpherePoints.end();)
	{
		Vec3 vPoint = it->first + vTargetEye;
		int& iType = it->second;

		Solution_t solution; CalculateAngle(tInfo.m_vLocalEye, vPoint, tInfo, iSimTime, solution, false);
		
		if (solution.m_iCalculated == CalculatedEnum::Bad)
			iType = 0;
		else if (abs(solution.m_flTime - TICKS_TO_TIME(iSimTime)) < tInfo.m_flRadiusTime || tInfo.m_flPrimeTime && iSimTime == tInfo.m_iPrimeTime)
			TracePoint(vPoint, iType, vTargetEye, tInfo, vPointDistances, checkPoint, i++);

		if (!(iType & ~PointTypeEnum::ObscuredMulti))
			it = vSpherePoints.erase(it);
		else
			++it;
	}

	std::sort(vPointDistances.begin(), vPointDistances.end(), [&](const auto& a, const auto& b) -> bool
		{
			return a.second < b.second;
		});

	std::vector<Point_t> vPoints = {};
	int iSplashCount = std::min(
		tInfo.m_flPrimeTime && iSimTime == tInfo.m_iPrimeTime ? Vars::Aimbot::Projectile::SplashCountDirect.Value : tInfo.m_iSplashCount,
		int(vPointDistances.size())
	);
	for (int i = 0; i < iSplashCount; i++)
		vPoints.push_back(vPointDistances[i].first);

	const Vec3 vOriginal = tTarget.m_pEntity->GetAbsOrigin();
	tTarget.m_pEntity->SetAbsOrigin(tTarget.m_vPos);
	for (auto it = vPoints.begin(); it != vPoints.end();)
	{
		auto& vPoint = *it;
		bool bValid = vPoint.m_Solution.m_iCalculated != CalculatedEnum::Pending;
		if (bValid)
		{
			Vec3 vPos; reinterpret_cast<CCollisionProperty*>(tTarget.m_pEntity->GetCollideable())->CalcNearestPoint(vPoint.m_vPoint, &vPos);
			bValid = vPoint.m_vPoint.DistTo(vPos) < tInfo.m_flRadius;
		}

		if (bValid)
			++it;
		else
			it = vPoints.erase(it);
	}
	tTarget.m_pEntity->SetAbsOrigin(vOriginal);

	return vPoints;
}

void CAimbotProjectile::SetupSplashPoints(Target_t& tTarget, std::vector<std::pair<Vec3, int>>& vSpherePoints, Info_t& tInfo, std::vector<std::pair<Vec3, Vec3>>& vSimplePoints)
{
	vSimplePoints.clear();
	Vec3 vTargetEye = tTarget.m_vPos + tInfo.m_vTargetEye;

	auto checkPoint = [&](CGameTrace& trace, bool& bErase, bool& bNormal)
		{
			bErase = !trace.m_pEnt || trace.fraction == 1.f || trace.surface.flags & SURF_SKY || !trace.m_pEnt->GetAbsVelocity().IsZero();
			if (bErase)
				return false;

			Point_t tPoint = { trace.endpos, {} };
			if (!tInfo.m_flGravity)
			{
				Vec3 vForward = (trace.endpos - tInfo.m_vLocalEye).Normalized();
				bNormal = vForward.Dot(trace.plane.normal) >= 0;
			}
			if (!bNormal)
			{
				CalculateAngle(tInfo.m_vLocalEye, tPoint.m_vPoint, tInfo, 0, tPoint.m_Solution, false);
				if (tInfo.m_flGravity)
				{
					Vec3 vPos = tInfo.m_vLocalEye + Vec3(0, 0, (tInfo.m_flGravity * 800.f * pow(tPoint.m_Solution.m_flTime, 2)) / 2);
					Vec3 vForward = (tPoint.m_vPoint - vPos).Normalized();
					bNormal = vForward.Dot(trace.plane.normal) >= 0;
				}
			}
			if (bNormal)
				return false;

			if (tPoint.m_Solution.m_iCalculated != CalculatedEnum::Bad)
			{
				vSimplePoints.emplace_back(tPoint.m_vPoint, trace.plane.normal);
				return true;
			}
			return false;
		};

	int i = 0;
	for (auto& vSpherePoint : vSpherePoints)
	{
		Vec3 vPoint = vSpherePoint.first + vTargetEye;
		int& iType = vSpherePoint.second;

		Solution_t solution; CalculateAngle(tInfo.m_vLocalEye, vPoint, tInfo, 0, solution, false);

		if (solution.m_iCalculated != CalculatedEnum::Bad)
			TracePoint(vPoint, iType, vTargetEye, tInfo, vSimplePoints, checkPoint, i++);
	}
}

std::vector<Point_t> CAimbotProjectile::GetSplashPointsSimple(Target_t& tTarget, std::vector<std::pair<Vec3, Vec3>>& vSpherePoints, Info_t& tInfo, int iSimTime)
{
	std::vector<std::pair<Point_t, float>> vPointDistances = {};

	Vec3 vTargetEye = tTarget.m_vPos + tInfo.m_vTargetEye;

	auto checkPoint = [&](Vec3& vPoint, bool& bErase)
		{
			Point_t tPoint = { vPoint, {} };
			CalculateAngle(tInfo.m_vLocalEye, tPoint.m_vPoint, tInfo, iSimTime, tPoint.m_Solution);

			bErase = tPoint.m_Solution.m_iCalculated == CalculatedEnum::Good;
			if (!bErase || !tInfo.m_flPrimeTime && int(tPoint.m_Solution.m_flTime / TICK_INTERVAL) + 1 != iSimTime)
				return false;

			vPointDistances.emplace_back(tPoint, tPoint.m_vPoint.DistTo(tTarget.m_vPos));
			return true;
		};
	for (auto it = vSpherePoints.begin(); it != vSpherePoints.end();)
	{
		Vec3& vPoint = it->first;
		Vec3& vNormal = it->second;
		bool bErase = false;

		if (checkPoint(vPoint, bErase))
		{
#ifdef SPLASH_DEBUG3
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bErase ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bErase ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(vNormal, vAngles);
				G::BoxStorage.emplace_back(vPoint, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::Halloween.Value, Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(tInfo.m_vLocalEye, vPoint), I::GlobalVars->curtime + 60.f, Vars::Colors::Halloween.Value);
			}
#endif

			/*
			CGameTrace trace = {};
			CTraceFilterWorldAndPropsOnly filter = {};

			SDK::Trace(vPoint + vNormal, vTargetEye, MASK_SHOT, &filter, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash trace check"]++;
#endif
#ifdef SPLASH_DEBUG3
			{
				Vec3 vMins = Vec3(-1, -1, -1) * (bErase ? 0.5f : 1.f), vMaxs = Vec3(1, 1, 1) * (bErase ? 0.5f : 1.f);
				Vec3 vAngles; Math::VectorAngles(trace.plane.normal, vAngles);
				G::BoxStorage.emplace_back(trace.endpos, vMins, vMaxs, vAngles, I::GlobalVars->curtime + 60.f, Vars::Colors::Powerup.Value.Alpha(Vars::Colors::Powerup.Value.a / (bErase ? 10 : 1)), Color_t(0, 0, 0, 0));
				G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(trace.startpos, trace.endpos), I::GlobalVars->curtime + 60.f, Vars::Colors::Powerup.Value.Alpha(Vars::Colors::Powerup.Value.a / (bErase ? 10 : 1)));
			}
#endif
			if (trace.fraction < 1.f)
				vPointDistances.pop_back();
			*/
		}

		if (bErase)
			it = vSpherePoints.erase(it);
		else
			++it;
	}

	std::sort(vPointDistances.begin(), vPointDistances.end(), [&](const auto& a, const auto& b) -> bool
		{
			return a.second < b.second;
		});

	std::vector<Point_t> vPoints = {};
	int iSplashCount = std::min(
		tInfo.m_flPrimeTime && iSimTime == tInfo.m_iPrimeTime ? Vars::Aimbot::Projectile::SplashCountDirect.Value : tInfo.m_iSplashCount,
		int(vPointDistances.size())
	);
	for (int i = 0; i < iSplashCount; i++)
		vPoints.push_back(vPointDistances[i].first);

	const Vec3 vOriginal = tTarget.m_pEntity->GetAbsOrigin();
	tTarget.m_pEntity->SetAbsOrigin(tTarget.m_vPos);
	for (auto it = vPoints.begin(); it != vPoints.end();)
	{
		auto& vPoint = *it;
		bool bValid = vPoint.m_Solution.m_iCalculated != CalculatedEnum::Pending;
		if (bValid)
		{
			Vec3 vPos = {}; reinterpret_cast<CCollisionProperty*>(tTarget.m_pEntity->GetCollideable())->CalcNearestPoint(vPoint.m_vPoint, &vPos);
			bValid = vPoint.m_vPoint.DistTo(vPos) < tInfo.m_flRadius;
		}

		if (bValid)
			++it;
		else
			it = vPoints.erase(it);
	}
	tTarget.m_pEntity->SetAbsOrigin(vOriginal);

	return vPoints;
}

static inline float AABBLine(Vec3 vMins, Vec3 vMaxs, Vec3 vStart, Vec3 vDir) noexcept
{
	// Enhanced AABB-ray intersection with improved numerical stability
	constexpr float kEpsilon = 1e-8f;
	
	float tMin = 0.0f;
	float tMax = std::numeric_limits<float>::max();
	
	// Process each axis with robust division handling
	for (int i = 0; i < 3; ++i) {
		const float origin = (&vStart.x)[i];
		const float direction = (&vDir.x)[i];
		const float boxMin = (&vMins.x)[i];
		const float boxMax = (&vMaxs.x)[i];
		
		if (std::abs(direction) < kEpsilon) {
			// Ray is parallel to the slab - check if origin is within bounds
			if (origin < boxMin || origin > boxMax) {
				return -1.0f; // No intersection
			}
		} else {
			// Calculate intersection distances with improved precision
			const float invDir = 1.0f / direction;
			float t1 = (boxMin - origin) * invDir;
			float t2 = (boxMax - origin) * invDir;
			
			// Ensure t1 <= t2
			if (t1 > t2) {
				std::swap(t1, t2);
			}
			
			// Update intersection interval
			tMin = std::max(tMin, t1);
			tMax = std::min(tMax, t2);
			
			// Early exit if no intersection possible
			if (tMin > tMax) {
				return -1.0f;
			}
		}
	}
	
	// Return the near intersection point (entry point)
	return std::max(0.0f, tMin);
}
static inline Vec3 PullPoint(Vec3 vPoint, Vec3 vLocalPos, Info_t& tInfo, Vec3 vMins, Vec3 vMaxs, Vec3 vTargetPos)
{
	auto HeightenLocalPos = [&]() -> Vec3
		{	// Enhanced trajectory calculation with improved numerical precision and projectile origin correction
			const double dGrav = static_cast<double>(tInfo.m_flGravity) * 800.0;
			if (dGrav < 1e-6)
				return vPoint;

			// Calculate preliminary angle to get weapon offset direction
			Vec3 vPrelimAngle = Math::CalcAngle(vLocalPos, vTargetPos);
			Vec3 vForward, vRight, vUp;
			Math::AngleVectors(vPrelimAngle, &vForward, &vRight, &vUp);
			
			// Calculate actual projectile spawn position with weapon offset
			Vec3 vProjectileOrigin = vLocalPos + (vForward * tInfo.m_vOffset.x) + (vRight * tInfo.m_vOffset.y) + (vUp * tInfo.m_vOffset.z);

			const Vec3 vDelta = vTargetPos - vProjectileOrigin;
			const double dDist = static_cast<double>(vDelta.Length2D());
			
			// Early exit for degenerate cases
			if (dDist < 1e-6)
				return vPoint;

			const double dVelocity = static_cast<double>(tInfo.m_flVelocity);
			const double dVelSq = dVelocity * dVelocity;
			const double dVelQuad = dVelSq * dVelSq;
			const double dDeltaZ = static_cast<double>(vDelta.z);
			
			// Enhanced discriminant calculation
			const double dRoot = dVelQuad - dGrav * (dGrav * dDist * dDist + 2.0 * dDeltaZ * dVelSq);
			if (dRoot < 0.0)
				return vPoint;
				
			const double dSqrtRoot = std::sqrt(dRoot);
			const double dPitch = std::atan((dVelSq - dSqrtRoot) / (dGrav * dDist));
			const double dCosPitch = std::cos(dPitch);
			
			if (std::abs(dCosPitch) < 1e-8)
				return vPoint;
				
			const double dTime = dDist / (dCosPitch * dVelocity) - static_cast<double>(tInfo.m_flOffsetTime);
			
			// Enhanced height calculation with improved precision
			const double dHeightOffset = (dGrav * dTime * dTime) * 0.5;
			return vProjectileOrigin + Vec3(0, 0, static_cast<float>(dHeightOffset));
		};

	Vec3 vProjectilePos = HeightenLocalPos();
	Vec3 vForward, vRight, vUp;
	Math::AngleVectors(Math::CalcAngle(vProjectilePos, vPoint), &vForward, &vRight, &vUp);
	
	// Enhanced AABB intersection calculation using corrected projectile position
	const Vec3 vDirection = vPoint - vProjectilePos;
	const float flDirectionLength = vDirection.Length();
	
	if (flDirectionLength < 1e-6f)
		return vProjectilePos;
		
	const Vec3 vNormalizedDirection = vDirection * (1.0f / flDirectionLength);
	const float flIntersectionT = AABBLine(vMins + vTargetPos, vMaxs + vTargetPos, vProjectilePos, vNormalizedDirection);
	
	// Clamp the intersection parameter to ensure valid results
	const float flClampedT = std::max(0.0f, std::min(flIntersectionT, flDirectionLength));
	
	return vProjectilePos + vNormalizedDirection * flClampedT;
}



static inline void SolveProjectileSpeed(CTFWeaponBase* pWeapon, const Vec3& vLocalPos, const Vec3& vTargetPos, float& flVelocity, float& flDragTime, const float flGravity) noexcept
{
	if (!F::ProjSim.obj->IsDragEnabled() || F::ProjSim.obj->m_dragBasis.IsZero()) [[likely]]
		return;

	// Enhanced numerical precision for drag calculations with improved distance accuracy
	const long double ldGrav = static_cast<long double>(flGravity) * 800.0L;
	const Vec3 vDelta = vTargetPos - vLocalPos;
	const long double ldDist = static_cast<long double>(vDelta.Length());  // Use 3D distance for better accuracy
	const long double ldVelocity = static_cast<long double>(flVelocity);
	const long double ldVelSq = ldVelocity * ldVelocity;

	// Early exit for degenerate cases with tighter tolerance
	if (ldDist < 1e-9L || ldVelocity < 1e-9L) [[unlikely]]
		return;

	// Enhanced ballistic calculation with improved numerical stability for demoman pipes
	const long double ldVelQuad = ldVelSq * ldVelSq;
	const long double ldDeltaZ = static_cast<long double>(vDelta.z);
	const long double ldDist2D = static_cast<long double>(vDelta.Length2D());
	
	// Use more accurate discriminant calculation for projectile physics
	const long double ldRoot = ldVelQuad - ldGrav * (ldGrav * ldDist2D * ldDist2D + 2.0L * ldDeltaZ * ldVelSq);
	
	if (ldRoot < 0.0L) [[unlikely]]
		return;

	const long double ldSqrtRoot = std::sqrt(ldRoot);
	const long double ldPitch = std::atan((ldVelSq - ldSqrtRoot) / (ldGrav * ldDist2D));
	const long double ldCosPitch = std::cos(ldPitch);
	
	if (std::abs(ldCosPitch) < 1e-12L) [[unlikely]]
		return;
		
	// More accurate time calculation using actual projectile path distance
	const long double ldTime = ldDist2D / (ldCosPitch * ldVelocity);
	const float flTime = static_cast<float>(ldTime);

	float flDrag = 0.0f;
	if (const float dragOverride = Vars::Aimbot::Projectile::DragOverride.Value; dragOverride > 0.0f) [[unlikely]]
	{
		flDrag = dragOverride;
	}
	else
	{
		// Get the weapon's item definition index
		const int weaponIndex = pWeapon->m_iItemDefinitionIndex();
		
		// Fixed drag calculation for demoman weapons - most should have zero drag for accurate trajectory
		struct DragEntry { int itemIndex; float minVel, maxVel, minDrag, maxDrag; };
		static constexpr std::array<DragEntry, 19> dragTable = {{
			{Demoman_m_GrenadeLauncher, 1200.0f, 1200.0f, 0.0f, 0.0f},  // Standard grenade launcher - no drag
			{Demoman_m_GrenadeLauncherR, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_FestiveGrenadeLauncher, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_Autumn, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_MacabreWeb, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_Rainbow, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_SweetDreams, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_CoffinNail, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_TopShelf, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_Warhawk, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_ButcherBird, 1200.0f, 1200.0f, 0.0f, 0.0f},
			{Demoman_m_TheIronBomber, 1200.0f, 1200.0f, 0.0f, 0.0f},  // Fixed: Iron Bomber should have no drag
			{Demoman_m_TheLochnLoad, 1513.0f, 1513.0f, 0.0f, 0.0f},   // Fixed: Loch-n-Load should have no drag
			{Demoman_m_TheLooseCannon, 1454.0f, 1454.0f, 0.0f, 0.0f}, // Fixed: Loose Cannon should have no drag
			{Demoman_s_StickybombLauncher, 900.0f, 2400.0f, 0.0f, 0.0f},  // Stickybomb launcher - no drag
			{Demoman_s_StickybombLauncherR, 900.0f, 2400.0f, 0.0f, 0.0f},
			{Demoman_s_FestiveStickybombLauncher, 900.0f, 2400.0f, 0.0f, 0.0f},
			{Demoman_s_TheQuickiebombLauncher, 900.0f, 2400.0f, 0.0f, 0.0f},
			{Demoman_s_TheScottishResistance, 900.0f, 2400.0f, 0.0f, 0.0f}  // Fixed: Scottish Resistance should have no drag
		}};
		
		// Fixed drag values for scout and other weapons - using constexpr for compile-time optimization
		struct FixedDragEntry { int itemIndex; float drag; };
		static constexpr std::array<FixedDragEntry, 9> fixedDragTable = {{
			{Scout_s_TheFlyingGuillotine, 0.0f},
			{Scout_s_TheFlyingGuillotineG, 0.310f},
			{Scout_t_TheSandman, 0.180f},
			{Scout_t_TheWrapAssassin, 0.285f},
			{Scout_s_MadMilk, 0.0f},
			{Scout_s_MutatedMilk, 0.0f},
			{Sniper_s_Jarate, 0.0f},
			{Sniper_s_FestiveJarate, 0.0f},
			{Sniper_s_TheSelfAwareBeautyMark, 0.057f}
		}};

		// Use std::find_if with lambda for better performance and readability
		const auto dragEntry = std::find_if(dragTable.begin(), dragTable.end(),
			[weaponIndex](const auto& entry) noexcept { return entry.itemIndex == weaponIndex; });
		
		if (dragEntry != dragTable.end() && dragEntry->maxVel > 0.0f) [[unlikely]]
		{
			// For weapons with fixed drag values, use the fixed value
			if (dragEntry->minVel == dragEntry->maxVel) {
				flDrag = dragEntry->minDrag;
			} else {
				flDrag = Math::RemapVal(flVelocity, dragEntry->minVel, dragEntry->maxVel, dragEntry->minDrag, dragEntry->maxDrag);
			}
		}
		else
		{
			// Check fixed drag table
			const auto fixedEntry = std::find_if(fixedDragTable.begin(), fixedDragTable.end(),
				[weaponIndex](const auto& entry) noexcept { return entry.itemIndex == weaponIndex; });
			
			if (fixedEntry != fixedDragTable.end()) [[unlikely]]
				flDrag = fixedEntry->drag;
		}
	}

	const float flOverride = Vars::Aimbot::Projectile::TimeOverride.Value;
	const float flDragDivisor = (flOverride > 0.0f) ? flOverride : 1.5f;
	
	// Enhanced drag integration using Runge-Kutta 4th order method for better accuracy
	constexpr int kDragSteps = 12; // Increased for better accuracy
	const long double ldStepTime = ldTime / static_cast<long double>(kDragSteps);
	long double ldIntegratedDrag = 0.0L;
	
	// Use Runge-Kutta 4th order for more accurate integration
	for (int i = 0; i < kDragSteps; ++i) {
		const long double ldCurrentTime = static_cast<long double>(i) * ldStepTime;
		const long double ldNextTime = static_cast<long double>(i + 1) * ldStepTime;
		const long double ldMidTime1 = ldCurrentTime + ldStepTime * 0.5L;
		const long double ldMidTime2 = ldCurrentTime + ldStepTime * 0.5L;
		
		// RK4 coefficients
		const long double ldDragCoeff = static_cast<long double>(flDrag);
		const long double k1 = ldVelocity * (1.0L - ldDragCoeff * ldCurrentTime) * ldDragCoeff;
		const long double k2 = ldVelocity * (1.0L - ldDragCoeff * ldMidTime1) * ldDragCoeff;
		const long double k3 = ldVelocity * (1.0L - ldDragCoeff * ldMidTime2) * ldDragCoeff;
		const long double k4 = ldVelocity * (1.0L - ldDragCoeff * ldNextTime) * ldDragCoeff;
		
		ldIntegratedDrag += (k1 + 2.0L * k2 + 2.0L * k3 + k4) * ldStepTime / 6.0L;
	}
	
	flDragTime = static_cast<float>(ldIntegratedDrag / static_cast<long double>(flDragDivisor));
	
	// Enhanced velocity calculation with exponential decay model for better accuracy
	const long double ldDragCoeff = static_cast<long double>(flDrag) * ldTime;
	flVelocity = static_cast<float>(ldVelocity * std::exp(-ldDragCoeff));
}
void CAimbotProjectile::CalculateAngle(const Vec3& vLocalPos, const Vec3& vTargetPos, Info_t& tInfo, int iSimTime, Solution_t& out, bool bAccuracy)
{
	if (out.m_iCalculated != CalculatedEnum::Pending)
		return;

	// Use long double for maximum precision in ballistic calculations
	const long double ldGrav = static_cast<long double>(tInfo.m_flGravity) * 800.0L;

	float flPitch, flYaw;
	{	// Enhanced trajectory calculation with improved numerical precision and projectile origin correction
		float flVelocity = tInfo.m_flVelocity, flDragTime = 0.f;
		
		// Calculate preliminary angle to get weapon offset direction
		Vec3 vPrelimAngle = Math::CalcAngle(vLocalPos, vTargetPos);
		Vec3 vForward, vRight, vUp;
		Math::AngleVectors(vPrelimAngle, &vForward, &vRight, &vUp);
		
		// Calculate actual projectile spawn position with corrected weapon offset
		// Fix coordinate transformation: offset.y should be negated for proper right vector
		Vec3 vProjectileOrigin = vLocalPos + (vForward * tInfo.m_vOffset.x) + (vRight * -tInfo.m_vOffset.y) + (vUp * tInfo.m_vOffset.z);
		
		// Use projectile origin for all calculations to fix alignment at long distances
		if (F::ProjSim.obj->IsDragEnabled() && !F::ProjSim.obj->m_dragBasis.IsZero())
		{
			SolveProjectileSpeed(tInfo.m_pWeapon, vProjectileOrigin, vTargetPos, flVelocity, flDragTime, tInfo.m_flGravity);
		}

		// Calculate delta from actual projectile spawn position, not eye position
		const Vec3 vDelta = vTargetPos - vProjectileOrigin;
		const long double ldDist2D = static_cast<long double>(vDelta.Length2D());
		const long double ldDist3D = static_cast<long double>(vDelta.Length());
		
		// Early exit for zero distance to avoid division by zero
		if (ldDist2D < 1e-9L) {
			out.m_iCalculated = CalculatedEnum::Bad;
			return;
		}

		// Calculate angle from projectile origin to target for proper alignment
		const Vec3 vAngleTo = Math::CalcAngle(vProjectileOrigin, vTargetPos);
		if (ldGrav < 1e-9L) // No gravity case
		{
			flPitch = -DEG2RAD(vAngleTo.x);
		}
		else
		{	// Enhanced ballistic trajectory calculation accounting for initial upward velocity
			const long double ldVelocity = static_cast<long double>(flVelocity);
			const long double ldVelSq = ldVelocity * ldVelocity;
			const long double ldDeltaZ = static_cast<long double>(vDelta.z);
			
			// Account for initial upward velocity for demoman projectiles (pipes/stickies get +200 units/sec upward)
			long double ldInitialUpVelocity = 0.0L;
			if (tInfo.m_pWeapon) {
				const int weaponID = tInfo.m_pWeapon->GetWeaponID();
				if (weaponID == TF_WEAPON_GRENADELAUNCHER || weaponID == TF_WEAPON_PIPEBOMBLAUNCHER || weaponID == TF_WEAPON_CANNON) {
					ldInitialUpVelocity = 200.0L; // Match ProjectileSimulation.cpp line 519/524
				}
			}
			
			// Enhanced ballistic formula accounting for initial upward velocity
			// Standard formula: tan(θ) = (v²±√(v⁴-g(gx²+2yv²))) / (gx)
			// Modified for initial upward velocity: account for vy₀ in height calculation
			const long double ldAdjustedDeltaZ = ldDeltaZ - (ldInitialUpVelocity * ldInitialUpVelocity) / (2.0L * ldGrav);
			const long double ldVelQuad = ldVelSq * ldVelSq;
			const long double ldDiscriminant = ldVelQuad - ldGrav * (ldGrav * ldDist2D * ldDist2D + 2.0L * ldAdjustedDeltaZ * ldVelSq);
			
			if (ldDiscriminant < 0.0L) {
				out.m_iCalculated = CalculatedEnum::Bad;
				return;
			}
			
			// Use the lower trajectory angle for better accuracy (more stable solution)
			const long double ldSqrtDiscriminant = std::sqrt(ldDiscriminant);
			const long double ldNumerator = ldVelSq - ldSqrtDiscriminant;
			const long double ldDenominator = ldGrav * ldDist2D;
			
			// Enhanced angle calculation with numerical stability check
			if (std::abs(ldDenominator) < 1e-15L) {
				out.m_iCalculated = CalculatedEnum::Bad;
				return;
			}
			
			const long double ldPitchRad = std::atan(ldNumerator / ldDenominator);
			flPitch = static_cast<float>(ldPitchRad);
			
			// Adjust pitch to account for initial upward velocity component
			if (ldInitialUpVelocity > 0.0L) {
				const long double ldUpwardAngle = std::atan(ldInitialUpVelocity / (ldVelocity * std::cos(ldPitchRad)));
				flPitch += static_cast<float>(ldUpwardAngle);
			}
			
			// Validate the calculated pitch angle
			if (!std::isfinite(flPitch) || std::abs(flPitch) > PI/2.0f) {
				out.m_iCalculated = CalculatedEnum::Bad;
				return;
			}
		}
		
		// Enhanced time calculation with improved precision using actual projectile path
		const long double ldCosPitch = std::cos(static_cast<long double>(flPitch));
		if (std::abs(ldCosPitch) < 1e-12L) {
			out.m_iCalculated = CalculatedEnum::Bad;
			return;
		}
		
		// Calculate flight time using 2D distance and pitch angle for accuracy
		const long double ldFlightTime = ldDist2D / (ldCosPitch * static_cast<long double>(flVelocity));
		out.m_flTime = static_cast<float>(ldFlightTime) - tInfo.m_flOffsetTime + flDragTime;
		out.m_flPitch = flPitch = -RAD2DEG(flPitch) - tInfo.m_vAngFix.x;
		out.m_flYaw = flYaw = vAngleTo.y - tInfo.m_vAngFix.y;
	}

	int iTimeTo = int(out.m_flTime / TICK_INTERVAL) + 1;
	if (out.m_iCalculated = iTimeTo > iSimTime ? CalculatedEnum::Time : CalculatedEnum::Pending)
		return;

	int iFlags = (bAccuracy ? ProjSimEnum::Trace : ProjSimEnum::None) | ProjSimEnum::NoRandomAngles | ProjSimEnum::PredictCmdNum;
#ifdef SPLASH_DEBUG6
	if (iFlags & ProjSimEnum::Trace)
	{
		if (Vars::Visuals::Trajectory::Override.Value)
		{
			if (!Vars::Visuals::Trajectory::Pipes.Value)
				mTraceCount["Setup trace calculate"]++;
		}
		else
		{
			switch (tInfo.m_pWeapon->GetWeaponID())
			{
			case TF_WEAPON_ROCKETLAUNCHER:
			case TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT:
			case TF_WEAPON_PARTICLE_CANNON:
			case TF_WEAPON_RAYGUN:
			case TF_WEAPON_DRG_POMSON:
			case TF_WEAPON_FLAREGUN:
			case TF_WEAPON_FLAREGUN_REVENGE:
			case TF_WEAPON_COMPOUND_BOW:
			case TF_WEAPON_CROSSBOW:
			case TF_WEAPON_SHOTGUN_BUILDING_RESCUE:
			case TF_WEAPON_SYRINGEGUN_MEDIC:
			case TF_WEAPON_FLAME_BALL:
				mTraceCount["Setup trace calculate"]++;
			}
		}
	}
#endif
	ProjectileInfo tProjInfo = {};
	if (out.m_iCalculated = !F::ProjSim.GetInfo(tInfo.m_pLocal, tInfo.m_pWeapon, { flPitch, flYaw, 0 }, tProjInfo, iFlags) ? CalculatedEnum::Bad : CalculatedEnum::Pending)
		return;

	{	// calculate trajectory from projectile origin with enhanced precision and proper alignment
		float flVelocity = tInfo.m_flVelocity, flDragTime = 0.f;
		SolveProjectileSpeed(tInfo.m_pWeapon, tProjInfo.m_vPos, vTargetPos, flVelocity, flDragTime, tInfo.m_flGravity);

		// Use the actual projectile spawn position for accurate calculations
		Vec3 vDelta = vTargetPos - tProjInfo.m_vPos;
		const long double ldDist2D = static_cast<long double>(vDelta.Length2D());
		const long double ldDist3D = static_cast<long double>(vDelta.Length());

		// Calculate angle from actual projectile position to target
		Vec3 vAngleTo = Math::CalcAngle(tProjInfo.m_vPos, vTargetPos);
		if (ldGrav < 1e-9L)
			out.m_flPitch = -DEG2RAD(vAngleTo.x);
		else
		{	// Enhanced ballistic calculation with improved precision for projectile origin
			const long double ldVelSq = static_cast<long double>(flVelocity * flVelocity);
			const long double ldVelQuad = ldVelSq * ldVelSq;
			const long double ldDeltaZ = static_cast<long double>(vDelta.z);
			
			const long double ldRoot = ldVelQuad - ldGrav * (ldGrav * ldDist2D * ldDist2D + 2.0L * ldDeltaZ * ldVelSq);
			if (ldRoot < 0.0L) {
				out.m_iCalculated = CalculatedEnum::Bad;
				return;
			}
			
			const long double ldNumerator = ldVelSq - std::sqrt(ldRoot);
			const long double ldDenominator = ldGrav * ldDist2D;
			
			if (std::abs(ldDenominator) < 1e-15L) {
				out.m_iCalculated = CalculatedEnum::Bad;
				return;
			}
			
			out.m_flPitch = static_cast<float>(std::atan(ldNumerator / ldDenominator));
		}
		
		// Enhanced time calculation with better precision
		const long double ldCosPitch = std::cos(static_cast<long double>(out.m_flPitch));
		if (std::abs(ldCosPitch) < 1e-12L) {
			out.m_iCalculated = CalculatedEnum::Bad;
			return;
		}
		
		out.m_flTime = static_cast<float>(ldDist2D / (ldCosPitch * static_cast<long double>(flVelocity))) + flDragTime;
	}

	{	// correct yaw with improved precision using proper projectile origin
		Vec3 vShootPos = (tProjInfo.m_vPos - vLocalPos).To2D();
		Vec3 vTarget = vTargetPos - vLocalPos;
		Vec3 vForward; Math::AngleVectors(tProjInfo.m_vAng, &vForward); vForward.Normalize2D();
		
		// Use double precision for quadratic solution
		const double dB = 2.0 * (static_cast<double>(vShootPos.x) * static_cast<double>(vForward.x) +
								 static_cast<double>(vShootPos.y) * static_cast<double>(vForward.y));
		const double dC = static_cast<double>(vShootPos.Length2DSqr()) - static_cast<double>(vTarget.Length2DSqr());
		
		auto vSolutions = Math::SolveQuadratic(1.0, dB, dC);
		if (!vSolutions.empty())
		{
			vShootPos += vForward * static_cast<float>(vSolutions.front());
			out.m_flYaw = flYaw - (RAD2DEG(atan2(vShootPos.y, vShootPos.x)) - flYaw);
			flYaw = RAD2DEG(atan2(vShootPos.y, vShootPos.x));
		}
	}

	{	// correct pitch with improved precision accounting for projectile origin offset
		if (ldGrav > 1e-9L)
		{
			flPitch -= tProjInfo.m_vAng.x;
			out.m_flPitch = -RAD2DEG(out.m_flPitch) + flPitch - tInfo.m_vAngFix.x;
		}
		else
		{
			// Account for projectile origin offset in pitch correction for non-gravity projectiles
			Vec3 vShootPos = Math::RotatePoint(tProjInfo.m_vPos - vLocalPos, {}, { 0, -flYaw, 0 }); vShootPos.y = 0;
			Vec3 vTarget = Math::RotatePoint(vTargetPos - vLocalPos, {}, { 0, -flYaw, 0 });
			Vec3 vForward; Math::AngleVectors(tProjInfo.m_vAng - Vec3(0, flYaw, 0), &vForward); vForward.y = 0; vForward.Normalize();
			
			// Use double precision for quadratic solution
			const double dB = 2.0 * (static_cast<double>(vShootPos.x) * static_cast<double>(vForward.x) +
									 static_cast<double>(vShootPos.z) * static_cast<double>(vForward.z));
			const double dC = (static_cast<double>(vShootPos.x * vShootPos.x) + static_cast<double>(vShootPos.z * vShootPos.z)) -
							  (static_cast<double>(vTarget.x * vTarget.x) + static_cast<double>(vTarget.z * vTarget.z));
			
			auto vSolutions = Math::SolveQuadratic(1.0, dB, dC);
			if (!vSolutions.empty())
			{
				vShootPos += vForward * static_cast<float>(vSolutions.front());
				out.m_flPitch = flPitch - (RAD2DEG(atan2(-vShootPos.z, vShootPos.x)) - flPitch);
			}
		}
	}

	iTimeTo = int(out.m_flTime / TICK_INTERVAL) + 1;
	out.m_iCalculated = iTimeTo > iSimTime ? CalculatedEnum::Time : CalculatedEnum::Good;
}



class CTraceFilterProjectileNoPlayer : public ITraceFilter
{
public:
	bool ShouldHitEntity(IHandleEntity* pServerEntity, int nContentsMask) override;
	TraceType_t GetTraceType() const override;
	CBaseEntity* pSkip = nullptr;
};
bool CTraceFilterProjectileNoPlayer::ShouldHitEntity(IHandleEntity* pServerEntity, int nContentsMask)
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
	case ETFClassID::CObjectSentrygun:
	case ETFClassID::CObjectDispenser:
	case ETFClassID::CObjectTeleporter: return true;
	}

	return false;
}
TraceType_t CTraceFilterProjectileNoPlayer::GetTraceType() const
{
	return TRACE_EVERYTHING;
}

// Optimized TestAngle function with modern C++ practices, const correctness, and performance improvements
bool CAimbotProjectile::TestAngle(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, Target_t& tTarget, Vec3& vPoint, Vec3& vAngles, int iSimTime, bool bSplash, bool* pHitSolid , std::vector<Vec3>* pProjectilePath) noexcept
{
	constexpr int iFlags = ProjSimEnum::Trace | ProjSimEnum::InitCheck | ProjSimEnum::NoRandomAngles | ProjSimEnum::PredictCmdNum;
	
#ifdef SPLASH_DEBUG6
	if (Vars::Visuals::Trajectory::Override.Value) [[unlikely]]
	{
		if (!Vars::Visuals::Trajectory::Pipes.Value)
			mTraceCount["Setup trace test"]++;
	}
	else
	{
		// Use constexpr array for better performance
		static constexpr std::array<int, 11> debugWeapons = {
			TF_WEAPON_ROCKETLAUNCHER, TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT, TF_WEAPON_PARTICLE_CANNON,
			TF_WEAPON_RAYGUN, TF_WEAPON_DRG_POMSON, TF_WEAPON_FLAREGUN, TF_WEAPON_FLAREGUN_REVENGE,
			TF_WEAPON_COMPOUND_BOW, TF_WEAPON_CROSSBOW, TF_WEAPON_SHOTGUN_BUILDING_RESCUE,
			TF_WEAPON_SYRINGEGUN_MEDIC, TF_WEAPON_FLAME_BALL
		};
		
		const int weaponID = pWeapon->GetWeaponID();
		if (std::find(debugWeapons.begin(), debugWeapons.end(), weaponID) != debugWeapons.end()) [[unlikely]]
			mTraceCount["Setup trace test"]++;
	}
	mTraceCount["Trace init check test"]++;
#endif

	ProjectileInfo tProjInfo = {};
	if (!F::ProjSim.GetInfo(pLocal, pWeapon, vAngles, tProjInfo, iFlags) || !F::ProjSim.Initialize(tProjInfo)) [[unlikely]]
		return false;

	bool bDidHit = false;
	CGameTrace trace = {};
	CTraceFilterProjectile filter = {};
	filter.pSkip = pLocal;
	CTraceFilterProjectileNoPlayer filterSplash = {};

#ifdef SPLASH_DEBUG5
	G::BoxStorage.emplace_back(vPoint, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, Vec3(), I::GlobalVars->curtime + 5.f, Color_t(255, 0, 0), Color_t(0, 0, 0, 0));
#endif

	// Early exit for non-gravity projectiles
	if (!tProjInfo.m_flGravity) [[unlikely]]
	{
		CTraceFilterWorldAndPropsOnly filterWorld = {};
		SDK::TraceHull(tProjInfo.m_vPos, vPoint, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, MASK_SOLID, &filterWorld, &trace);
#ifdef SPLASH_DEBUG6
		mTraceCount["Nograv trace"]++;
#endif
#ifdef SPLASH_DEBUG5
		G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(tProjInfo.m_vPos, vPoint), I::GlobalVars->curtime + 5.f, Color_t(0, 0, 0));
#endif
		if (trace.fraction < 0.999f) [[unlikely]]
			return false;
	}

	// Cache original position and set target position
	const Vec3 vOriginal = tTarget.m_pEntity->GetAbsOrigin();
	tTarget.m_pEntity->SetAbsOrigin(tTarget.m_vPos);
	
	// Pre-calculate constants for better performance
	const int weaponID = pWeapon->GetWeaponID();
	const bool bIsLunchbox = (weaponID == TF_WEAPON_LUNCHBOX);
	const int splashInterval = Vars::Aimbot::Projectile::SplashTraceInterval.Value;
	
	bool bPrimeTime = false;
	Vec3 vStaticPos = {};
	
	// Main simulation loop with optimized branching
	for (int n = 1; n <= iSimTime; ++n)
	{
		const Vec3 vOld = F::ProjSim.GetOrigin();
		F::ProjSim.RunTick(tProjInfo);
		const Vec3 vNew = F::ProjSim.GetOrigin();

		if (bDidHit) [[unlikely]]
		{
			trace.endpos = vNew;
			continue;
		}

		// Optimized tracing logic with branch prediction hints
		if (!bSplash) [[likely]]
		{
			SDK::TraceHull(vOld, vNew, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, MASK_SOLID, &filter, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Direct trace"]++;
#endif
#ifdef SPLASH_DEBUG5
			G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(vOld, vNew), I::GlobalVars->curtime + 5.f, Color_t(255, 0, 0));
#endif
		}
		else
		{
			if (n == 1 || bPrimeTime) [[unlikely]]
				vStaticPos = vOld;
			if (n % splashInterval && n != iSimTime && !bPrimeTime) [[likely]]
				continue;

			SDK::TraceHull(vStaticPos, vNew, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, MASK_SOLID, &filterSplash, &trace);
#ifdef SPLASH_DEBUG6
			mTraceCount["Splash trace"]++;
#endif
#ifdef SPLASH_DEBUG5
			G::LineStorage.emplace_back(std::pair<Vec3, Vec3>(vStaticPos, vNew), I::GlobalVars->curtime + 5.f, Color_t(255, 0, 0));
#endif
			vStaticPos = vNew;
		}
		
		if (trace.DidHit()) [[unlikely]]
		{
			if (pHitSolid) [[likely]]
				*pHitSolid = true;

			// Optimized validation logic with branch prediction
			const bool bTime = bSplash ?
				trace.endpos.DistTo(vPoint) < tProjInfo.m_flVelocity * TICK_INTERVAL + tProjInfo.m_vHull.z :
				iSimTime - n < 5 || bIsLunchbox;
			const bool bTarget = trace.m_pEnt == tTarget.m_pEntity || bSplash;
			bool bValid = bTarget && bTime;
			
			if (bValid && bSplash) [[unlikely]]
			{
				bValid = SDK::VisPosWorld(nullptr, tTarget.m_pEntity, trace.endpos, vPoint, MASK_SOLID);
#ifdef SPLASH_DEBUG6
				mTraceCount["Splash vispos"]++;
#endif
				if (bValid) [[likely]]
				{
					// Use constexpr array for rocket weapons check
					static constexpr std::array<int, 3> rocketWeapons = {
						TF_WEAPON_ROCKETLAUNCHER, TF_WEAPON_ROCKETLAUNCHER_DIRECTHIT, TF_WEAPON_PARTICLE_CANNON
					};
					
					if (std::find(rocketWeapons.begin(), rocketWeapons.end(), weaponID) != rocketWeapons.end()) [[unlikely]]
					{
						CGameTrace eyeTrace = {};
						CTraceFilterWorldAndPropsOnly filter = {};
						SDK::Trace(trace.endpos + trace.plane.normal, tTarget.m_vPos + tTarget.m_pEntity->As<CTFPlayer>()->GetViewOffset(), MASK_SHOT, &filter, &eyeTrace);
						bValid = eyeTrace.fraction == 1.f;
#ifdef SPLASH_DEBUG6
						mTraceCount["Rocket trace"]++;
#endif
					}
				}
			}

#ifdef SPLASH_DEBUG5
			G::BoxStorage.pop_back();
			if (bValid)
				G::BoxStorage.emplace_back(vPoint, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, Vec3(), I::GlobalVars->curtime + 5.f, Color_t(0, 255, 0), Color_t(0, 0, 0, 0));
			else if (!bTime)
			{
				G::BoxStorage.emplace_back(vPoint, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, Vec3(), I::GlobalVars->curtime + 5.f, Color_t(255, 0, 255), Color_t(0, 0, 0, 0));
				if (bSplash)
				{
					G::BoxStorage.emplace_back(trace.endpos, Vec3(-1, -1, -1), Vec3(1, 1, 1), Vec3(), I::GlobalVars->curtime + 5.f, Color_t(0, 0, 0), Color_t(0, 0, 0, 0));
					G::BoxStorage.emplace_back(vPoint, Vec3(-1, -1, -1), Vec3(1, 1, 1), Vec3(), I::GlobalVars->curtime + 5.f, Color_t(255, 255, 255), Color_t(0, 0, 0, 0));
				}
			}
			else
				G::BoxStorage.emplace_back(vPoint, tProjInfo.m_vHull * -1, tProjInfo.m_vHull, Vec3(), I::GlobalVars->curtime + 5.f, Color_t(0, 0, 255), Color_t(0, 0, 0, 0));
#endif

			if (bValid) [[likely]]
			{
				if (bSplash) [[unlikely]]
				{
					const int iPopCount = splashInterval - static_cast<int>(trace.fraction * splashInterval);
					for (int i = 0; i < iPopCount && !tProjInfo.m_vPath.empty(); ++i)
						tProjInfo.m_vPath.pop_back();
				}

				// Optimized aim type checking
				const int aimType = Vars::Aimbot::General::AimType.Value;
				if ((aimType == Vars::Aimbot::General::AimTypeEnum::Smooth ||
					 aimType == Vars::Aimbot::General::AimTypeEnum::Assistive) &&
					tTarget.m_nAimedHitbox == HITBOX_HEAD) [[unlikely]]
				{
					const Vec3 vOffset = (trace.endpos - vNew) + (vOriginal - tTarget.m_vPos);

					Vec3 vOld = F::ProjSim.GetOrigin() + vOffset;
					F::ProjSim.RunTick(tProjInfo);
					Vec3 vNew = F::ProjSim.GetOrigin() + vOffset;

					CGameTrace boneTrace = {};
					SDK::Trace(vOld, vNew, MASK_SHOT, &filter, &boneTrace);
#ifdef SPLASH_DEBUG6
					mTraceCount["Huntsman trace"]++;
#endif
					boneTrace.endpos -= vOffset;

					if (boneTrace.DidHit() && (boneTrace.m_pEnt != tTarget.m_pEntity || boneTrace.hitbox != HITBOX_HEAD))
						goto continueLoop;

					if (!boneTrace.DidHit()) // loop and see if closest hitbox is head
					{
						const auto pModel = tTarget.m_pEntity->GetModel();
						if (!pModel) goto continueLoop;
						const auto pHDR = I::ModelInfoClient->GetStudiomodel(pModel);
						if (!pHDR) goto continueLoop;
						const auto pSet = pHDR->pHitboxSet(tTarget.m_pEntity->As<CTFPlayer>()->m_nHitboxSet());
						if (!pSet) goto continueLoop;

						const auto pBones = H::Entities.GetBones(tTarget.m_pEntity->entindex());
						if (!pBones) goto continueLoop;

						Vec3 vForward = vOld - vNew;
						vForward.Normalize();
						const Vec3 vPos = boneTrace.endpos + vForward * 16.0f + vOriginal - tTarget.m_vPos;

						float closestDist = std::numeric_limits<float>::max();
						int closestId = -1;
						
						for (int i = 0; i < pSet->numhitboxes; ++i)
						{
							const auto pBox = pSet->pHitbox(i);
							if (!pBox) continue;

							Vec3 vCenter;
							Math::VectorTransform((pBox->bbmin + pBox->bbmax) * 0.5f, pBones[pBox->bone], vCenter);

							const float flDist = vPos.DistTo(vCenter);
							if (flDist < closestDist)
							{
								closestDist = flDist;
								closestId = i;
							}
						}

						if (closestId != 0) goto continueLoop;
						bDidHit = true;
					}
				}

				bDidHit = true;
			}
			else if (!bSplash && bTarget && weaponID == TF_WEAPON_PIPEBOMBLAUNCHER) [[unlikely]]
			{
				// run for more ticks to check for splash
				iSimTime = n + 5;
				bSplash = bPrimeTime = true;
			}
			else
				break;

			continueLoop:
			if (!bSplash) [[likely]]
				trace.endpos = vNew;

			if (!bTarget || (bSplash && !bPrimeTime)) [[unlikely]]
				break;
		}
	}
	
	// Restore original position
	tTarget.m_pEntity->SetAbsOrigin(vOriginal);

	if (bDidHit && pProjectilePath) [[likely]]
	{
		tProjInfo.m_vPath.push_back(trace.endpos);
		*pProjectilePath = std::move(tProjInfo.m_vPath);
	}

	return bDidHit;
}



int CAimbotProjectile::CanHit(Target_t& tTarget, CTFPlayer* pLocal, CTFWeaponBase* pWeapon,
	std::vector<Vec3>* pPlayerPath, std::vector<Vec3>* pProjectilePath, std::vector<DrawBox_t>* pBoxes, float* pTimeTo)
{
	//if (Vars::Aimbot::General::Ignore.Value & Vars::Aimbot::General::IgnoreEnum::Unsimulated && H::Entities.GetChoke(tTarget.m_pEntity->entindex()) > Vars::Aimbot::General::TickTolerance.Value)
	//	return false;

	PlayerStorage tStorage;
	ProjectileInfo tProjInfo = {};

	int iMaxTime, iSplash; Info_t tInfo = { pLocal, pWeapon };
	const float flSize = tTarget.m_pEntity->GetSize().Length();
	{
		int iFlags = ProjSimEnum::NoRandomAngles | ProjSimEnum::PredictCmdNum;
		if (!F::ProjSim.GetInfo(pLocal, pWeapon, {}, tProjInfo, iFlags) || !F::ProjSim.Initialize(tProjInfo, false))
			return false;

		tInfo.m_vLocalEye = pLocal->GetShootPos();
		tInfo.m_vTargetEye = tTarget.m_pEntity->As<CTFPlayer>()->GetViewOffset();
		F::MoveSim.Initialize(tTarget.m_pEntity, tStorage);
		tTarget.m_vPos = tTarget.m_pEntity->m_vecOrigin();

		tInfo.m_flLatency = F::Backtrack.GetReal() + TICKS_TO_TIME(F::Backtrack.GetAnticipatedChoke());

		iMaxTime = TIME_TO_TICKS(std::min(tProjInfo.m_flLifetime, Vars::Aimbot::Projectile::MaxSimulationTime.Value));

		Vec3 vVelocity = F::ProjSim.GetVelocity();
		
		// For demoman projectiles, use horizontal velocity for ballistic calculations
		// since the upward component is handled separately in the physics simulation
		if (pWeapon) {
			const int weaponID = pWeapon->GetWeaponID();
			if (weaponID == TF_WEAPON_GRENADELAUNCHER || weaponID == TF_WEAPON_PIPEBOMBLAUNCHER || weaponID == TF_WEAPON_CANNON) {
				// Use 2D velocity magnitude for ballistic calculations to match the physics simulation
				tInfo.m_flVelocity = vVelocity.Length2D();
			} else {
				tInfo.m_flVelocity = vVelocity.Length();
			}
		} else {
			tInfo.m_flVelocity = vVelocity.Length();
		}
		
		Math::VectorAngles(vVelocity, tInfo.m_vAngFix);

		tInfo.m_vHull = tProjInfo.m_vHull.Min(3);
		tInfo.m_vOffset = tProjInfo.m_vPos - tInfo.m_vLocalEye; tInfo.m_vOffset.y *= -1;
		tInfo.m_flOffsetTime = tInfo.m_vOffset.Length() / tInfo.m_flVelocity; // silly

		tInfo.m_flGravity = tProjInfo.m_flGravity;
		tInfo.m_flRadius = GetSplashRadius(pWeapon, pLocal); tInfo.m_flRadiusTime = tInfo.m_flRadius / tInfo.m_flVelocity;
		tInfo.m_flBoundingTime = tInfo.m_flRadiusTime + flSize / tInfo.m_flVelocity;

		iSplash = Vars::Aimbot::Projectile::SplashPrediction.Value && tInfo.m_flRadius ? Vars::Aimbot::Projectile::SplashPrediction.Value : Vars::Aimbot::Projectile::SplashPredictionEnum::Off;
		tInfo.m_iSplashCount = !tProjInfo.m_flGravity ? Vars::Aimbot::Projectile::SplashCountDirect.Value : Vars::Aimbot::Projectile::SplashCountArc.Value;

		tInfo.m_iSplashMode = GetSplashMode(pWeapon);
		tInfo.m_flPrimeTime = PrimeTime(pWeapon);
		tInfo.m_iPrimeTime = TIME_TO_TICKS(tInfo.m_flPrimeTime);
	}

	int iReturn = false;

	int iMulti = Vars::Aimbot::Projectile::SplashMode.Value;

	auto mDirectPoints = iSplash == Vars::Aimbot::Projectile::SplashPredictionEnum::Only ? std::unordered_map<int, Vec3>() : GetDirectPoints(tTarget, tInfo);
	auto vSpherePoints = !iSplash ? std::vector<std::pair<Vec3, int>>() : ComputeSphere(tInfo.m_flRadius + flSize, Vars::Aimbot::Projectile::SplashPoints.Value, Vars::Aimbot::Projectile::SplashNthRoot.Value);
#ifdef SPLASH_DEBUG4
	for (auto& [vPoint, _] : vSpherePoints)
		G::BoxStorage.emplace_back(tTarget.m_pEntity->m_vecOrigin() + tInfo.m_vTargetEye + vPoint * tInfo.m_flRadius / (tInfo.m_flRadius + flSize), Vec3(-1, -1, -1), Vec3(1, 1, 1), Vec3(), I::GlobalVars->curtime + 60.f, Color_t(0, 0, 0, 0), Vars::Colors::Local.Value);
#endif

	
	Vec3 vAngleTo, vPredicted, vTarget;
	int iLowestPriority = std::numeric_limits<int>::max(); float flLowestDist = std::numeric_limits<float>::max();
	int iLowestSmoothPriority = iLowestPriority; float flLowestSmoothDist = flLowestDist;
	
	const int iStartTime = 1 - TIME_TO_TICKS(tInfo.m_flLatency);
	const float flInvMaxTime = 1.f / static_cast<float>(iMaxTime - iStartTime);
	
	for (int i = iStartTime; i <= iMaxTime; i++)
	{
		if (!tStorage.m_bFailed)
		{
			F::MoveSim.RunTick(tStorage);
			tTarget.m_vPos = tStorage.m_vPredictedOrigin;
		}
		if (i < 0)
			continue;

		bool bDirectBreaks = true;
		std::vector<Point_t> vSplashPoints = {};
		if (iSplash)
		{
			Solution_t solution; CalculateAngle(tInfo.m_vLocalEye, tTarget.m_vPos, tInfo, i, solution, false);
			if (solution.m_iCalculated != CalculatedEnum::Bad)
			{
				bDirectBreaks = false;

				const float flTimeTo = solution.m_flTime - TICKS_TO_TIME(i);
				if (flTimeTo < tInfo.m_flBoundingTime)
				{
					static std::vector<std::pair<Vec3, Vec3>> vSimplePoints = {};
					if (iMulti == Vars::Aimbot::Projectile::SplashModeEnum::Single)
					{
						SetupSplashPoints(tTarget, vSpherePoints, tInfo, vSimplePoints);
						if (!vSimplePoints.empty())
							iMulti++;
						else
						{
							iSplash = Vars::Aimbot::Projectile::SplashPredictionEnum::Off;
							goto skipSplash;
						}
					}

					if ((iMulti == Vars::Aimbot::Projectile::SplashModeEnum::Multi ? vSpherePoints.empty() : vSimplePoints.empty())
						|| flTimeTo < -tInfo.m_flBoundingTime && (tInfo.m_flPrimeTime ? i > tInfo.m_iPrimeTime : true))
						break;
					else if (tInfo.m_flPrimeTime ? i >= tInfo.m_iPrimeTime : true)
					{
						if (iMulti == Vars::Aimbot::Projectile::SplashModeEnum::Multi)
							vSplashPoints = GetSplashPoints(tTarget, vSpherePoints, tInfo, i);
						else
							vSplashPoints = GetSplashPointsSimple(tTarget, vSimplePoints, tInfo, i);
					}
				}
			}
		}
		skipSplash:
		if (bDirectBreaks && mDirectPoints.empty())
			break;

		std::vector<std::tuple<Point_t, int, int>> vPoints = {};
		for (auto& [iIndex, vPoint] : mDirectPoints)
		{
			Point_t tPoint = { tTarget.m_vPos + vPoint, {} };
			vPoints.emplace_back(tPoint, iIndex + (iSplash == Vars::Aimbot::Projectile::SplashPredictionEnum::Prefer ? tInfo.m_iSplashCount : 0), iIndex);
		}
		for (auto& vPoint : vSplashPoints)
			vPoints.emplace_back(vPoint, iSplash == Vars::Aimbot::Projectile::SplashPredictionEnum::Include ? 3 : 0, -1);

		int j = 0;
		for (auto& [vPoint, iPriority, iIndex] : vPoints) // get most ideal point
		{
			const bool bSplash = iIndex == -1;
			Vec3 vOriginalPoint = vPoint.m_vPoint;

			if (Vars::Aimbot::Projectile::HuntsmanPullPoint.Value && tTarget.m_nAimedHitbox == HITBOX_HEAD)
				vPoint.m_vPoint = PullPoint(vPoint.m_vPoint, tInfo.m_vLocalEye, tInfo, tTarget.m_pEntity->m_vecMins() + tProjInfo.m_vHull, tTarget.m_pEntity->m_vecMaxs() - tProjInfo.m_vHull, tTarget.m_vPos);
				//vPoint.m_vPoint = PullPoint(vPoint.m_vPoint, tInfo.m_vLocalEye, tInfo, tTarget.m_pEntity->m_vecMins(), tTarget.m_pEntity->m_vecMaxs(), tTarget.m_vPos);

			float flDist = bSplash ? tTarget.m_vPos.DistTo(vPoint.m_vPoint) : flLowestDist;
			bool bPriority = bSplash ? iPriority <= iLowestPriority : iPriority < iLowestPriority;
			bool bTime = bSplash || tInfo.m_iPrimeTime < i || tStorage.m_MoveData.m_vecVelocity.IsZero();
			bool bDist = !bSplash || flDist < flLowestDist;
			if (!bSplash && !bPriority)
				mDirectPoints.erase(iIndex);
			if (!bPriority || !bTime || !bDist)
				continue;

			CalculateAngle(tInfo.m_vLocalEye, vPoint.m_vPoint, tInfo, i, vPoint.m_Solution);
			if (!bSplash && (vPoint.m_Solution.m_iCalculated == CalculatedEnum::Good || vPoint.m_Solution.m_iCalculated == CalculatedEnum::Bad))
				mDirectPoints.erase(iIndex);
			if (vPoint.m_Solution.m_iCalculated != CalculatedEnum::Good)
				continue;

			if (Vars::Aimbot::Projectile::HuntsmanPullPoint.Value && tTarget.m_nAimedHitbox == HITBOX_HEAD)
			{
				Solution_t tSolution;
				CalculateAngle(tInfo.m_vLocalEye, vOriginalPoint, tInfo, std::numeric_limits<int>::max(), tSolution);
				vPoint.m_Solution.m_flPitch = tSolution.m_flPitch, vPoint.m_Solution.m_flYaw = tSolution.m_flYaw;
			}

			Vec3 vAngles; Aim(G::CurrentUserCmd->viewangles, { vPoint.m_Solution.m_flPitch, vPoint.m_Solution.m_flYaw, 0.f }, vAngles);
			std::vector<Vec3> vProjLines; bool bHitSolid = false;

			if (TestAngle(pLocal, pWeapon, tTarget, vPoint.m_vPoint, vAngles, i, bSplash, &bHitSolid, &vProjLines))
			{
				iLowestPriority = iPriority; flLowestDist = flDist;
				vAngleTo = vAngles, vPredicted = tTarget.m_vPos, vTarget = vOriginalPoint;
				*pTimeTo = vPoint.m_Solution.m_flTime + tInfo.m_flLatency;
				*pPlayerPath = tStorage.m_vPath;
				pPlayerPath->push_back(tStorage.m_MoveData.m_vecAbsOrigin);
				*pProjectilePath = vProjLines;
			}
			else switch (Vars::Aimbot::General::AimType.Value)
			{
			case Vars::Aimbot::General::AimTypeEnum::Smooth:
			case Vars::Aimbot::General::AimTypeEnum::Assistive:
			{
				if (/*Vars::Aimbot::General::AssistStrength.Value != 0.f &&*/ Vars::Aimbot::General::AssistStrength.Value != 100.f)
				{
					bPriority = bSplash ? iPriority <= iLowestSmoothPriority : iPriority < flLowestSmoothDist;
					bDist = !bSplash || flDist < flLowestDist;
					if (!bPriority || !bDist)
						continue;

					Vec3 vPlainAngles; Aim({}, { vPoint.m_Solution.m_flPitch, vPoint.m_Solution.m_flYaw, 0.f }, vPlainAngles, Vars::Aimbot::General::AimTypeEnum::Plain);
					if (TestAngle(pLocal, pWeapon, tTarget, vPoint.m_vPoint, vPlainAngles, i, bSplash, &bHitSolid))
					{
						iLowestSmoothPriority = iPriority; flLowestSmoothDist = flDist;
						vAngleTo = vAngles, vPredicted = tTarget.m_vPos;
						*pPlayerPath = tStorage.m_vPath;
						pPlayerPath->push_back(tStorage.m_MoveData.m_vecAbsOrigin);
						iReturn = 2;
					}
				}
			}
			}

			if (!j && bHitSolid)
				*pTimeTo = vPoint.m_Solution.m_flTime + tInfo.m_flLatency;
			j++;
		}
	}
	F::MoveSim.Restore(tStorage);

	tTarget.m_vPos = vTarget;

	if (tTarget.m_iTargetType != TargetEnum::Player || !tStorage.m_bFailed) // don't attempt to aim at players when movesim fails
	{
		bool bMain = iLowestPriority != std::numeric_limits<int>::max();
		bool bAny = bMain || iLowestSmoothPriority != std::numeric_limits<int>::max();

		tTarget.m_vAngleTo = vAngleTo;

		if (bAny && (Vars::Colors::BoundHitboxEdge.Value.a || Vars::Colors::BoundHitboxFace.Value.a || Vars::Colors::BoundHitboxEdgeClipped.Value.a || Vars::Colors::BoundHitboxFaceClipped.Value.a))
		{
			tInfo.m_vHull = tInfo.m_vHull.Max(1);
			float flProjectileTime = TICKS_TO_TIME(pProjectilePath->size());
			float flTargetTime = tStorage.m_bFailed ? flProjectileTime : TICKS_TO_TIME(pPlayerPath->size());

			if (Vars::Visuals::Hitbox::BoundsEnabled.Value & Vars::Visuals::Hitbox::BoundsEnabledEnum::OnShot)
			{
				if (Vars::Colors::BoundHitboxEdge.Value.a || Vars::Colors::BoundHitboxFace.Value.a)
					pBoxes->emplace_back(vPredicted, tTarget.m_pEntity->m_vecMins(), tTarget.m_pEntity->m_vecMaxs(), Vec3(), I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flTargetTime : Vars::Visuals::Hitbox::DrawDuration.Value), Vars::Colors::BoundHitboxEdge.Value, Vars::Colors::BoundHitboxFace.Value);
				if (Vars::Colors::BoundHitboxEdgeClipped.Value.a || Vars::Colors::BoundHitboxFaceClipped.Value.a)
					pBoxes->emplace_back(vPredicted, tTarget.m_pEntity->m_vecMins(), tTarget.m_pEntity->m_vecMaxs(), Vec3(), I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flTargetTime : Vars::Visuals::Hitbox::DrawDuration.Value), Vars::Colors::BoundHitboxEdgeClipped.Value, Vars::Colors::BoundHitboxFaceClipped.Value, true);
			}

			if (bMain && Vars::Visuals::Hitbox::BoundsEnabled.Value & Vars::Visuals::Hitbox::BoundsEnabledEnum::AimPoint)
			{
				if (Vars::Colors::BoundHitboxEdge.Value.a || Vars::Colors::BoundHitboxFace.Value.a)
					pBoxes->emplace_back(vTarget, tInfo.m_vHull * -1, tInfo.m_vHull, Vec3(), I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flProjectileTime : Vars::Visuals::Hitbox::DrawDuration.Value), Vars::Colors::BoundHitboxEdge.Value, Vars::Colors::BoundHitboxFace.Value);
				if (Vars::Colors::BoundHitboxEdgeClipped.Value.a || Vars::Colors::BoundHitboxFaceClipped.Value.a)
					pBoxes->emplace_back(vTarget, tInfo.m_vHull * -1, tInfo.m_vHull, Vec3(), I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flProjectileTime : Vars::Visuals::Hitbox::DrawDuration.Value), Vars::Colors::BoundHitboxEdgeClipped.Value, Vars::Colors::BoundHitboxFaceClipped.Value, true);

				/*
				if (Vars::Debug::Info.Value && tTarget.m_nAimedHitbox == HITBOX_HEAD) // huntsman head
				{
					const Vec3 vOriginOffset = tTarget.m_pEntity->m_vecOrigin() - vPredicted;

					auto pBones = H::Entities.GetBones(tTarget.m_pEntity->entindex());
					if (!pBones)
						return true;

					auto vBoxes = F::Visuals.GetHitboxes(pBones, tTarget.m_pEntity->As<CTFPlayer>(), { HITBOX_HEAD });
					for (auto& bBox : vBoxes)
					{
						bBox.m_vPos -= vOriginOffset;
						bBox.m_flTime = I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flTargetTime : Vars::Visuals::Hitbox::DrawDuration.Value);
						pBoxes->push_back(bBox);
					}
				}
				if (Vars::Debug::Info.Value && tTarget.m_nAimedHitbox == HITBOX_HEAD) // huntsman head, broken; removeme once 254 is fixed
				{
					const Vec3 vOriginOffset = tTarget.m_pEntity->m_vecOrigin() - vPredicted;

					auto pBones = H::Entities.GetBones(tTarget.m_pEntity->entindex());
					if (!pBones)
						return true;

					auto vBoxes = F::Visuals.GetHitboxes(pBones, tTarget.m_pEntity->As<CTFPlayer>(), { HITBOX_HEAD });
					for (auto& bBox : vBoxes)
					{
						bBox.m_vPos -= vOriginOffset;
						bBox.m_flTime = I::GlobalVars->curtime + (Vars::Visuals::Simulation::Timed.Value ? flTargetTime : Vars::Visuals::Hitbox::DrawDuration.Value);
						bBox.m_vRotation = Vec3();
						pBoxes->push_back(bBox);
					}
				}
				*/
			}
		}

		if (bMain)
			return true;
	}

	return iReturn;
}



bool CAimbotProjectile::Aim(Vec3 vCurAngle, Vec3 vToAngle, Vec3& vOut, int iMethod)
{
	if (Vec3* pDoubletapAngle = F::Ticks.GetShootAngle())
	{
		vOut = *pDoubletapAngle;
		return true;
	}

	Math::ClampAngles(vToAngle);

	switch (iMethod)
	{
	case Vars::Aimbot::General::AimTypeEnum::Plain:
	case Vars::Aimbot::General::AimTypeEnum::Silent:
	case Vars::Aimbot::General::AimTypeEnum::Locking:
		vOut = vToAngle;
		return false;
	case Vars::Aimbot::General::AimTypeEnum::Smooth:
		vOut = vCurAngle.LerpAngle(vToAngle, Vars::Aimbot::General::AssistStrength.Value / 100.f);
		return true;
	case Vars::Aimbot::General::AimTypeEnum::Assistive:
		Vec3 vMouseDelta = G::CurrentUserCmd->viewangles.DeltaAngle(G::LastUserCmd->viewangles);
		Vec3 vTargetDelta = vToAngle.DeltaAngle(G::LastUserCmd->viewangles);
		float flMouseDelta = vMouseDelta.Length2D(), flTargetDelta = vTargetDelta.Length2D();
		vTargetDelta = vTargetDelta.Normalized() * std::min(flMouseDelta, flTargetDelta);
		vOut = vCurAngle - vMouseDelta + vMouseDelta.LerpAngle(vTargetDelta, Vars::Aimbot::General::AssistStrength.Value / 100.f);
		return true;
	}

	return false;
}

// assume angle calculated outside with other overload
void CAimbotProjectile::Aim(CUserCmd* pCmd, Vec3& vAngle)
{
	switch (Vars::Aimbot::General::AimType.Value)
	{
	case Vars::Aimbot::General::AimTypeEnum::Plain:
	case Vars::Aimbot::General::AimTypeEnum::Smooth:
	case Vars::Aimbot::General::AimTypeEnum::Assistive:
		//pCmd->viewangles = vAngle; // retarded, overshooting with this uncommented
		I::EngineClient->SetViewAngles(vAngle);
		break;
	case Vars::Aimbot::General::AimTypeEnum::Silent:
	{
		bool bDoubleTap = F::Ticks.m_bDoubletap || F::Ticks.GetTicks(H::Entities.GetWeapon()) || F::Ticks.m_bSpeedhack;
		auto pWeapon = H::Entities.GetWeapon();
		if (G::Attacking == 1 || bDoubleTap || pWeapon && pWeapon->GetWeaponID() == TF_WEAPON_FLAMETHROWER)
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

static inline void CancelShot(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd, int& iLastTickCancel)
{
	switch (pWeapon->GetWeaponID())
	{
	case TF_WEAPON_COMPOUND_BOW:
	{
		pCmd->buttons |= IN_ATTACK2;
		pCmd->buttons &= ~IN_ATTACK;
		break;
	}
	case TF_WEAPON_CANNON:
	case TF_WEAPON_PIPEBOMBLAUNCHER:
	{
		for (int i = 0; i < MAX_WEAPONS; i++)
		{
			auto pSwap = pLocal->GetWeaponFromSlot(i);
			if (!pSwap || pSwap == pWeapon || !pSwap->CanBeSelected())
				continue;

			pCmd->weaponselect = pSwap->entindex();
			iLastTickCancel = pWeapon->entindex();
			break;
		}
	}
	}
}

bool CAimbotProjectile::RunMain(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	const int nWeaponID = pWeapon->GetWeaponID();

	static int iStaticAimType = Vars::Aimbot::General::AimType.Value;
	const int iLastAimType = iStaticAimType;
	const int iRealAimType = Vars::Aimbot::General::AimType.Value;

	switch (nWeaponID)
	{
	case TF_WEAPON_COMPOUND_BOW:
	case TF_WEAPON_PIPEBOMBLAUNCHER:
	case TF_WEAPON_CANNON:
		if (!Vars::Aimbot::General::AutoShoot.Value && G::Attacking && !iRealAimType && iLastAimType)
			Vars::Aimbot::General::AimType.Value = iLastAimType;
		break;
	default:
		if (G::Throwing && !iRealAimType && iLastAimType)
			Vars::Aimbot::General::AimType.Value = iLastAimType;
	}
	iStaticAimType = Vars::Aimbot::General::AimType.Value;

	if (F::AimbotGlobal.ShouldHoldAttack(pWeapon))
		pCmd->buttons |= IN_ATTACK;
	if (!Vars::Aimbot::General::AimType.Value
		|| !F::AimbotGlobal.ShouldAim() && nWeaponID != TF_WEAPON_PIPEBOMBLAUNCHER && nWeaponID != TF_WEAPON_CANNON && nWeaponID != TF_WEAPON_FLAMETHROWER)
		return false;

	auto vTargets = SortTargets(pLocal, pWeapon);
	if (vTargets.empty())
		return false;

	if (Vars::Aimbot::Projectile::Modifiers.Value & Vars::Aimbot::Projectile::ModifiersEnum::ChargeWeapon && iRealAimType
		&& (nWeaponID == TF_WEAPON_COMPOUND_BOW || nWeaponID == TF_WEAPON_PIPEBOMBLAUNCHER))
	{
		pCmd->buttons |= IN_ATTACK;
		if (!G::CanPrimaryAttack && !G::Reloading && Vars::Aimbot::General::AimType.Value == Vars::Aimbot::General::AimTypeEnum::Silent)
			return false;
	}

	if (!G::AimTarget.m_iEntIndex)
		G::AimTarget = { vTargets.front().m_pEntity->entindex(), I::GlobalVars->tickcount, 0 };

#if defined(SPLASH_DEBUG1) || defined(SPLASH_DEBUG2) || defined(SPLASH_DEBUG3) || defined(SPLASH_DEBUG5)
	G::LineStorage.clear();
#endif
#if defined(SPLASH_DEBUG1) || defined(SPLASH_DEBUG2) || defined(SPLASH_DEBUG3) || defined(SPLASH_DEBUG4) || defined(SPLASH_DEBUG5)
	G::BoxStorage.clear();
#endif
	for (auto& tTarget : vTargets)
	{
		float flTimeTo = 0.f; std::vector<Vec3> vPlayerPath, vProjectilePath; std::vector<DrawBox_t> vBoxes = {};
		const int iResult = CanHit(tTarget, pLocal, pWeapon, &vPlayerPath, &vProjectilePath, &vBoxes, &flTimeTo);
		if (iResult != 1 && pWeapon->GetWeaponID() == TF_WEAPON_CANNON && Vars::Aimbot::Projectile::Modifiers.Value & Vars::Aimbot::Projectile::ModifiersEnum::ChargeWeapon && !(pCmd->buttons & IN_ATTACK))
		{
			float flCharge = pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() > 0.f
				? std::clamp(pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() - I::GlobalVars->curtime, 0.f, 1.f)
				: 1.f;
			if (!flTimeTo)
				flTimeTo = std::numeric_limits<float>::max();
			if (flCharge < flTimeTo)
			{
				if (pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() > 0.f)
					CancelShot(pLocal, pWeapon, pCmd, m_iLastTickCancel);
			}
			else
			{
				if (m_iLastTickCancel)
					pCmd->weaponselect = m_iLastTickCancel = 0;
				pCmd->buttons |= IN_ATTACK;
			}
		}
		if (!iResult) continue;
		if (iResult == 2)
		{
			G::AimTarget = { tTarget.m_pEntity->entindex(), I::GlobalVars->tickcount, 0 };

			bool bPlayerPath = Vars::Visuals::Simulation::PlayerPath.Value;
			bool bBoxes = Vars::Visuals::Hitbox::BoundsEnabled.Value & (Vars::Visuals::Hitbox::BoundsEnabledEnum::OnShot | Vars::Visuals::Hitbox::BoundsEnabledEnum::AimPoint);
			if (bPlayerPath || bBoxes)
			{
				G::PathStorage.clear();
				G::BoxStorage.clear();
				G::LineStorage.clear();

				if (bPlayerPath)
				{
					if (Vars::Colors::PlayerPath.Value.a)
						G::PathStorage.emplace_back(vPlayerPath, Vars::Visuals::Simulation::Timed.Value ? -int(vPlayerPath.size()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
					if (Vars::Colors::PlayerPathClipped.Value.a)
						G::PathStorage.emplace_back(vPlayerPath, Vars::Visuals::Simulation::Timed.Value ? -int(vPlayerPath.size()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
				}
				if (bBoxes)
					G::BoxStorage.insert(G::BoxStorage.end(), vBoxes.begin(), vBoxes.end());
			}

			Aim(pCmd, tTarget.m_vAngleTo);
			break;
		}

		G::AimTarget = { tTarget.m_pEntity->entindex(), I::GlobalVars->tickcount };
		G::AimPoint = { tTarget.m_vPos, I::GlobalVars->tickcount };

		if (Vars::Aimbot::General::AutoShoot.Value)
		{
			switch (nWeaponID)
			{
			case TF_WEAPON_COMPOUND_BOW:
			case TF_WEAPON_PIPEBOMBLAUNCHER:
				pCmd->buttons |= IN_ATTACK;
				if (pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() > 0.f)
					pCmd->buttons &= ~IN_ATTACK;
				break;
			case TF_WEAPON_CANNON:
				pCmd->buttons |= IN_ATTACK;
				if (pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() > 0.f)
				{
					if (m_iLastTickCancel)
						pCmd->weaponselect = m_iLastTickCancel = 0;
					if (Vars::Aimbot::Projectile::Modifiers.Value & Vars::Aimbot::Projectile::ModifiersEnum::ChargeWeapon)
					{
						float flCharge = pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() - I::GlobalVars->curtime;
						if (std::clamp(flCharge, 0.f, 1.f) < flTimeTo)
							pCmd->buttons &= ~IN_ATTACK;
					}
					else
						pCmd->buttons &= ~IN_ATTACK;
				}
				break;
			case TF_WEAPON_BAT_WOOD:
			case TF_WEAPON_BAT_GIFTWRAP:
			case TF_WEAPON_LUNCHBOX:
				pCmd->buttons &= ~IN_ATTACK, pCmd->buttons |= IN_ATTACK2;
				break;
			default:
				pCmd->buttons |= IN_ATTACK;
				if (pWeapon->m_iItemDefinitionIndex() == Soldier_m_TheBeggarsBazooka)
				{
					if (pWeapon->m_iClip1() > 0)
						pCmd->buttons &= ~IN_ATTACK;
				}
			}
		}

		F::Aimbot.m_bRan = G::Attacking = SDK::IsAttacking(pLocal, pWeapon, pCmd, true);

		if (G::Attacking == 1 || !Vars::Aimbot::General::AutoShoot.Value)
		{
			bool bPlayerPath = Vars::Visuals::Simulation::PlayerPath.Value;
			bool bProjectilePath = Vars::Visuals::Simulation::ProjectilePath.Value && (G::Attacking == 1 || Vars::Debug::Info.Value);
			bool bBoxes = Vars::Visuals::Hitbox::BoundsEnabled.Value & (Vars::Visuals::Hitbox::BoundsEnabledEnum::OnShot | Vars::Visuals::Hitbox::BoundsEnabledEnum::AimPoint);
			if (bPlayerPath || bProjectilePath || bBoxes)
			{
				G::PathStorage.clear();
				G::BoxStorage.clear();
				G::LineStorage.clear();

				if (bPlayerPath)
				{
					if (Vars::Colors::PlayerPath.Value.a)
						G::PathStorage.emplace_back(vPlayerPath, Vars::Visuals::Simulation::Timed.Value ? -int(vPlayerPath.size()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPath.Value, Vars::Visuals::Simulation::PlayerPath.Value);
					if (Vars::Colors::PlayerPathClipped.Value.a)
						G::PathStorage.emplace_back(vPlayerPath, Vars::Visuals::Simulation::Timed.Value ? -int(vPlayerPath.size()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::PlayerPathClipped.Value, Vars::Visuals::Simulation::PlayerPath.Value, true);
				}
				if (bProjectilePath)
				{
					if (Vars::Colors::ProjectilePath.Value.a)
						G::PathStorage.emplace_back(vProjectilePath, Vars::Visuals::Simulation::Timed.Value ? -int(vProjectilePath.size()) - TIME_TO_TICKS(F::Backtrack.GetReal()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::ProjectilePath.Value, Vars::Visuals::Simulation::ProjectilePath.Value);
					if (Vars::Colors::ProjectilePathClipped.Value.a)
						G::PathStorage.emplace_back(vProjectilePath, Vars::Visuals::Simulation::Timed.Value ? -int(vProjectilePath.size()) - TIME_TO_TICKS(F::Backtrack.GetReal()) : I::GlobalVars->curtime + Vars::Visuals::Simulation::DrawDuration.Value, Vars::Colors::ProjectilePathClipped.Value, Vars::Visuals::Simulation::ProjectilePath.Value, true);
				}
				if (bBoxes)
					G::BoxStorage.insert(G::BoxStorage.end(), vBoxes.begin(), vBoxes.end());
			}
		}

		Aim(pCmd, tTarget.m_vAngleTo);
		if (G::PSilentAngles)
		{
			switch (nWeaponID)
			{
			case TF_WEAPON_FLAMETHROWER: // angles show up anyways
			case TF_WEAPON_CLEAVER: // can't psilent with these weapons, they use SetContextThink
			case TF_WEAPON_JAR:
			case TF_WEAPON_JAR_MILK:
			case TF_WEAPON_JAR_GAS:
			case TF_WEAPON_BAT_WOOD:
			case TF_WEAPON_BAT_GIFTWRAP:
				G::PSilentAngles = false, G::SilentAngles = true;
			}
		}
		return true;
	}

	return false;
}

void CAimbotProjectile::Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	const bool bSuccess = RunMain(pLocal, pWeapon, pCmd);
#if defined(SPLASH_DEBUG6)
	if (Vars::Aimbot::General::AimType.Value && !mTraceCount.empty())
	{
		int iTraceCount = 0;
		for (auto& [_, iTraces] : mTraceCount)
			iTraceCount += iTraces;
		SDK::Output("Traces", std::format("{}", iTraceCount).c_str());
		for (auto& [sType, iTraces] : mTraceCount)
			SDK::Output("Traces", std::format("{}: {}", sType, iTraces).c_str());
	}
	mTraceCount.clear();
#endif

	float flAmount = 0.f;
	if (pWeapon->GetWeaponID() == TF_WEAPON_PIPEBOMBLAUNCHER)
	{
		const float flCharge = pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() > 0.f ? I::GlobalVars->curtime - pWeapon->As<CTFPipebombLauncher>()->m_flChargeBeginTime() : 0.f;
		flAmount = Math::RemapVal(flCharge, 0.f, SDK::AttribHookValue(4.f, "stickybomb_charge_rate", pWeapon), 0.f, 1.f);
	}
	else if (pWeapon->GetWeaponID() == TF_WEAPON_CANNON)
	{
		const float flMortar = SDK::AttribHookValue(0.f, "grenade_launcher_mortar_mode", pWeapon);
		const float flCharge = pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() > 0.f ? I::GlobalVars->curtime - pWeapon->As<CTFGrenadeLauncher>()->m_flDetonateTime() : -flMortar;
		flAmount = flMortar ? Math::RemapVal(flCharge, -flMortar, 0.f, 0.f, 1.f) : 0.f;
	}

	if (pWeapon->GetWeaponID() == TF_WEAPON_PIPEBOMBLAUNCHER && G::OriginalMove.m_iButtons & IN_ATTACK && Vars::Aimbot::Projectile::AutoRelease.Value && flAmount > Vars::Aimbot::Projectile::AutoRelease.Value / 100)
		pCmd->buttons &= ~IN_ATTACK;
	else if (G::CanPrimaryAttack && Vars::Aimbot::Projectile::Modifiers.Value & Vars::Aimbot::Projectile::ModifiersEnum::CancelCharge)
	{
		if (m_bLastTickHeld && (G::LastUserCmd->buttons & IN_ATTACK && !(pCmd->buttons & IN_ATTACK) && !bSuccess || flAmount > 0.95f))
			CancelShot(pLocal, pWeapon, pCmd, m_iLastTickCancel);
	}

	m_bLastTickHeld = Vars::Aimbot::General::AimType.Value;
}
