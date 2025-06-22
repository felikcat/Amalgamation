// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "AimbotGlobal.h"

#include "../Aimbot.h"
#include "../../Players/PlayerUtils.h"
#include "../../Ticks/Ticks.h"

void CAimbotGlobal::SortTargets(std::vector<Target_t>& vTargets, const int iMethod) noexcept
{	// Sort by preference
	std::sort(vTargets.begin(), vTargets.end(), [iMethod](const Target_t& a, const Target_t& b) noexcept -> bool
		{
			switch (iMethod)
			{
			case Vars::Aimbot::General::TargetSelectionEnum::FOV: return a.m_flFOVTo < b.m_flFOVTo;
			case Vars::Aimbot::General::TargetSelectionEnum::Distance: return a.m_flDistTo < b.m_flDistTo;
			default: return false;
			}
		});
}

void CAimbotGlobal::SortPriority(std::vector<Target_t>& vTargets) noexcept
{	// Sort by priority
	std::sort(vTargets.begin(), vTargets.end(), [](const Target_t& a, const Target_t& b) noexcept -> bool
		{
			return a.m_nPriority > b.m_nPriority;
		});
}

// this won't prevent shooting bones outside of fov
bool CAimbotGlobal::PlayerBoneInFOV(CTFPlayer* pTarget, const Vec3 vLocalPos, const Vec3 vLocalAngles, float& flFOVTo, Vec3& vPos, Vec3& vAngleTo, const int iHitboxes) noexcept
{
	matrix3x4 aBones[MAXSTUDIOBONES];
	if (!pTarget->SetupBones(aBones, MAXSTUDIOBONES, BONE_USED_BY_ANYTHING, I::GlobalVars->curtime)) 
		return false;

	float flMinFOV = 180.f;
	const int numHitboxes = pTarget->GetNumOfHitboxes();
	const uint32_t modelHash = H::Entities.GetModel(pTarget->entindex());
	const float aimFOV = Vars::Aimbot::General::AimFOV.Value;
	
	for (int nHitbox = 0; nHitbox < numHitboxes; ++nHitbox)
	{
		if (!IsHitboxValid(modelHash, nHitbox, iHitboxes)) 
			continue;

		const Vec3 vCurPos = pTarget->GetHitboxCenter(aBones, nHitbox);
		const Vec3 vCurAngleTo = Math::CalcAngle(vLocalPos, vCurPos);
		const float flCurFOVTo = Math::CalcFov(vLocalAngles, vCurAngleTo);

		if (flCurFOVTo < flMinFOV) 
		{
			vPos = vCurPos;
			vAngleTo = vCurAngleTo;
			flFOVTo = flMinFOV = flCurFOVTo;
		}
	}

	return flMinFOV < aimFOV;
}

bool CAimbotGlobal::IsHitboxValid(const uint32_t uHash, const int nHitbox, const int iHitboxes) noexcept
{
	// Early return for special case
	if (nHitbox == -1) 
		return true;

	// Use constexpr for compile-time hash
	constexpr uint32_t saxtonHaleHash = FNV1A::Hash32Const("models/vsh/player/saxton_hale.mdl");
	
	if (uHash == saxtonHaleHash) 
	{
		switch (nHitbox)
		{
		case HITBOX_SAXTON_HEAD:
			return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Head;
		case HITBOX_SAXTON_BODY:
		case HITBOX_SAXTON_THORAX:
		case HITBOX_SAXTON_CHEST:
		case HITBOX_SAXTON_UPPER_CHEST:
		case HITBOX_SAXTON_NECK:
			return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Body;
		case HITBOX_SAXTON_PELVIS:
			return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Pelvis;
		case HITBOX_SAXTON_LEFT_UPPER_ARM:
		case HITBOX_SAXTON_LEFT_FOREARM:
		case HITBOX_SAXTON_LEFT_HAND:
		case HITBOX_SAXTON_RIGHT_UPPER_ARM:
		case HITBOX_SAXTON_RIGHT_FOREARM:
		case HITBOX_SAXTON_RIGHT_HAND:
			return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Arms;
		case HITBOX_SAXTON_LEFT_THIGH:
		case HITBOX_SAXTON_LEFT_CALF:
		case HITBOX_SAXTON_LEFT_FOOT:
		case HITBOX_SAXTON_RIGHT_THIGH:
		case HITBOX_SAXTON_RIGHT_CALF:
		case HITBOX_SAXTON_RIGHT_FOOT:
			return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Legs;
		default:
			return false;
		}
	}
	
	// Standard hitbox handling
	switch (nHitbox)
	{
	case HITBOX_HEAD:
		return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Head;
	case HITBOX_BODY:
	case HITBOX_THORAX:
	case HITBOX_CHEST:
	case HITBOX_UPPER_CHEST:
		return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Body;
	case HITBOX_PELVIS:
		return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Pelvis;
	case HITBOX_LEFT_UPPER_ARM:
	case HITBOX_LEFT_FOREARM:
	case HITBOX_LEFT_HAND:
	case HITBOX_RIGHT_UPPER_ARM:
	case HITBOX_RIGHT_FOREARM:
	case HITBOX_RIGHT_HAND:
		return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Arms;
	case HITBOX_LEFT_THIGH:
	case HITBOX_LEFT_CALF:
	case HITBOX_LEFT_FOOT:
	case HITBOX_RIGHT_THIGH:
	case HITBOX_RIGHT_CALF:
	case HITBOX_RIGHT_FOOT:
		return iHitboxes & Vars::Aimbot::Hitscan::HitboxesEnum::Legs;
	default:
		return false;
	}
}

bool CAimbotGlobal::ShouldIgnore(CBaseEntity* pEntity, CTFPlayer* pLocal, CTFWeaponBase* pWeapon) noexcept
{
	if (pEntity->IsDormant()) 
		return true;

	if (const auto pGameRules = I::TFGameRules()) 
	{
		if (pGameRules->m_bTruceActive() && pLocal->m_iTeamNum() != pEntity->m_iTeamNum()) 
			return true;
	}

	const auto classID = pEntity->GetClassID();
	switch (classID)
	{
	case ETFClassID::CTFPlayer:
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		if (pPlayer == pLocal || !pPlayer->IsAlive() || pPlayer->IsAGhost()) 
			return true;

		if (pLocal->m_iTeamNum() == pEntity->m_iTeamNum()) 
			return false;

		if (F::PlayerUtils.IsIgnored(pPlayer->entindex())) 
			return true;

		// Cache ignore flags for better performance
		const int ignoreFlags = Vars::Aimbot::General::Ignore.Value;
		const int playerIndex = pPlayer->entindex();
		
		if ((ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Friends && H::Entities.IsFriend(playerIndex))
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Party && H::Entities.InParty(playerIndex))
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Invulnerable && pPlayer->IsInvulnerable() && SDK::AttribHookValue(0, "crit_forces_victim_to_laugh", pWeapon) <= 0)
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Cloaked && pPlayer->m_flInvisibility() && pPlayer->m_flInvisibility() >= Vars::Aimbot::General::IgnoreCloak.Value / 100.f)
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::DeadRinger && pPlayer->m_bFeignDeathReady())
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Taunting && pPlayer->IsTaunting())
			|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Disguised && pPlayer->InCond(TF_COND_DISGUISED))) 
			return true;
		if (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Vaccinator) 
		{
			switch (G::PrimaryWeaponType)
			{
			case EWeaponType::HITSCAN:
				if (pPlayer->InCond(TF_COND_MEDIGUN_UBER_BULLET_RESIST) && SDK::AttribHookValue(0, "mod_pierce_resists_absorbs", pWeapon) != 0) 
					return true;
				break;
			case EWeaponType::PROJECTILE:
				{
					const auto weaponID = pWeapon->GetWeaponID();
					switch (weaponID)
					{
					case TF_WEAPON_FLAMETHROWER:
					case TF_WEAPON_FLAREGUN:
						if (pPlayer->InCond(TF_COND_MEDIGUN_UBER_FIRE_RESIST)) 
							return true;
						break;
					case TF_WEAPON_COMPOUND_BOW:
						if (pPlayer->InCond(TF_COND_MEDIGUN_UBER_BULLET_RESIST)) 
							return true;
						break;
					default:
						if (pPlayer->InCond(TF_COND_MEDIGUN_UBER_BLAST_RESIST)) 
							return true;
					}
				}
				break;
			}
		}

		return false;
	}
	case ETFClassID::CObjectSentrygun:
	case ETFClassID::CObjectDispenser:
	case ETFClassID::CObjectTeleporter:
	{
		const auto pBuilding = pEntity->As<CBaseObject>();
		const int targetFlags = Vars::Aimbot::General::Target.Value;

		if ((!(targetFlags & Vars::Aimbot::General::TargetEnum::Sentry) && pBuilding->IsSentrygun())
			|| (!(targetFlags & Vars::Aimbot::General::TargetEnum::Dispenser) && pBuilding->IsDispenser())
			|| (!(targetFlags & Vars::Aimbot::General::TargetEnum::Teleporter) && pBuilding->IsTeleporter())) 
			return true;

		if (pLocal->m_iTeamNum() == pEntity->m_iTeamNum()) 
			return false;

		if (const auto pOwner = pBuilding->m_hBuilder().Get()) 
		{
			const int ownerIndex = pOwner->entindex();
			if (F::PlayerUtils.IsIgnored(ownerIndex)) 
				return true;

			const int ignoreFlags = Vars::Aimbot::General::Ignore.Value;
			if ((ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Friends && H::Entities.IsFriend(ownerIndex))
				|| (ignoreFlags & Vars::Aimbot::General::IgnoreEnum::Party && H::Entities.InParty(ownerIndex))) 
				return true;
		}

		return false;
	}
	case ETFClassID::CTFGrenadePipebombProjectile:
	{
		if (!(Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::Stickies)) 
			return true;

		if (pLocal->m_iTeamNum() == pEntity->m_iTeamNum()) 
			return true;

		const auto pProjectile = pEntity->As<CTFGrenadePipebombProjectile>();
		if (const auto pOwner = pProjectile->m_hThrower().Get()) 
		{
			if (F::PlayerUtils.IsIgnored(pOwner->entindex())) 
				return true;
		}

		if (pProjectile->m_iType() != TF_GL_MODE_REMOTE_DETONATE || !pProjectile->m_bTouched()) 
			return true;

		return false;
	}
	case ETFClassID::CEyeballBoss:
	case ETFClassID::CHeadlessHatman:
	case ETFClassID::CMerasmus:
	case ETFClassID::CTFBaseBoss:
	{
		if (!(Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::NPCs))
			return true;

		if (pEntity->m_iTeamNum() != TF_TEAM_HALLOWEEN)
			return true;

		return false;
	}
	case ETFClassID::CTFTankBoss:
	case ETFClassID::CZombie:
	{
		if (!(Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::NPCs))
			return true;

		if (pLocal->m_iTeamNum() == pEntity->m_iTeamNum())
			return true;

		return false;
	}
	case ETFClassID::CTFPumpkinBomb:
	case ETFClassID::CTFGenericBomb:
	{
		if (!(Vars::Aimbot::General::Target.Value & Vars::Aimbot::General::TargetEnum::Bombs))
			return true;

		return false;
	}
	}

	return true;
}

int CAimbotGlobal::GetPriority(const int iIndex) noexcept
{
	return F::PlayerUtils.GetPriority(iIndex);
}

bool CAimbotGlobal::ShouldAim() noexcept
{
	const auto aimType = Vars::Aimbot::General::AimType.Value;
	switch (aimType)
	{
	case Vars::Aimbot::General::AimTypeEnum::Plain:
	case Vars::Aimbot::General::AimTypeEnum::Silent:
		// for performance reasons, the F::Ticks.m_bDoubletap condition is not a great fix here
		// and actually properly predicting when shots will be fired should likely be done over this, but it's fine for now
		if (!G::CanPrimaryAttack && !G::Reloading && !F::Ticks.m_bDoubletap && !F::Ticks.m_bSpeedhack) 
			return false;
		break;
	default:
		break;
	}

	return true;
}

bool CAimbotGlobal::ShouldHoldAttack(CTFWeaponBase* pWeapon) noexcept
{
	const auto holdsFire = Vars::Aimbot::General::AimHoldsFire.Value;
	switch (holdsFire)
	{
	case Vars::Aimbot::General::AimHoldsFireEnum::MinigunOnly:
		if (pWeapon->GetWeaponID() != TF_WEAPON_MINIGUN) 
			break;
		[[fallthrough]];
	case Vars::Aimbot::General::AimHoldsFireEnum::Always:
		if (!F::Aimbot.m_bRunningSecondary && !G::CanPrimaryAttack && (G::LastUserCmd->buttons & IN_ATTACK) && Vars::Aimbot::General::AimType.Value && !pWeapon->IsInReload()) 
			return true;
		break;
	default:
		break;
	}
	return false;
}

// will not predict for projectile weapons
bool CAimbotGlobal::ValidBomb(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CBaseEntity* pBomb) noexcept
{
	if (G::PrimaryWeaponType == EWeaponType::PROJECTILE) 
		return false;

	const Vec3 vOrigin = pBomb->m_vecOrigin();
	constexpr float bombRadius = 300.f;
	const int targetFlags = Vars::Aimbot::General::Target.Value;
	const int localTeam = pLocal->m_iTeamNum();

	CBaseEntity* pEntity;
	for (CEntitySphereQuery sphere(vOrigin, bombRadius);
		(pEntity = sphere.GetCurrentEntity()) != nullptr;
		sphere.NextEntity())
	{
		if (!pEntity || pEntity == pLocal || pEntity->m_iTeamNum() == localTeam) 
			continue;

		// Early check for player validity
		if (pEntity->IsPlayer()) 
		{
			const auto pPlayer = pEntity->As<CTFPlayer>();
			if (!pPlayer->IsAlive() || pPlayer->IsAGhost()) 
				continue;
		}

		Vec3 vPos;
		reinterpret_cast<CCollisionProperty*>(pEntity->GetCollideable())->CalcNearestPoint(vOrigin, &vPos);
		if (vOrigin.DistTo(vPos) > bombRadius) 
			continue;

		// Cache target type checks
		const bool isPlayer = pEntity->IsPlayer() && (targetFlags & Vars::Aimbot::General::TargetEnum::Players);
		const bool isSentry = pEntity->IsSentrygun() && (targetFlags & Vars::Aimbot::General::TargetEnum::Sentry);
		const bool isDispenser = pEntity->IsDispenser() && (targetFlags & Vars::Aimbot::General::TargetEnum::Dispenser);
		const bool isTeleporter = pEntity->IsTeleporter() && (targetFlags & Vars::Aimbot::General::TargetEnum::Teleporter);
		const bool isNPC = pEntity->IsNPC() && (targetFlags & Vars::Aimbot::General::TargetEnum::NPCs);
		
		if (isPlayer || isSentry || isDispenser || isTeleporter || isNPC) 
		{
			if (isPlayer && ShouldIgnore(pEntity->As<CTFPlayer>(), pLocal, pWeapon)) 
				continue;

			const Vec3 targetPos = isPlayer ? pEntity->m_vecOrigin() + pEntity->As<CTFPlayer>()->GetViewOffset() : pEntity->GetCenter();
			if (!SDK::VisPosProjectile(pBomb, pEntity, vOrigin, targetPos, MASK_SHOT)) 
				continue;

			return true;
		}
	}

	return false;
}