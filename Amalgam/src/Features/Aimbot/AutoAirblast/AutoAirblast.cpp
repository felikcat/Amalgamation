// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "AutoAirblast.h"

#include "../../Backtrack/Backtrack.h"

static inline bool ShouldTargetProjectile(CBaseEntity* pProjectile, CTFPlayer* pLocal)
{
	if (pProjectile->m_iTeamNum() == pLocal->m_iTeamNum())
		return false;

	switch (pProjectile->GetClassID())
	{
	case ETFClassID::CTFGrenadePipebombProjectile:
		if (pProjectile->As<CTFGrenadePipebombProjectile>()->m_bTouched())
			return false;
	}

	return true;
}

static inline Vec3 PredictOrigin(Vec3& vOrigin, Vec3 vVelocity, float flLatency, bool bTrace = true, Vec3 vMins = {}, Vec3 vMaxs = {}, unsigned int nMask = MASK_SOLID)
{
	if (vVelocity.IsZero() || !flLatency)
		return vOrigin;

	Vec3 vTo = vOrigin + vVelocity * flLatency;
	if (!bTrace)
		return vTo;

	CGameTrace trace = {};
	CTraceFilterWorldAndPropsOnly filter = {};

	SDK::TraceHull(vOrigin, vTo, vMins, vMaxs, nMask, &filter, &trace);
	return vOrigin + (vTo - vOrigin) * trace.fraction;
}

bool CAutoAirblast::CanAirblastEntity(CTFPlayer* pLocal, CBaseEntity* pEntity, Vec3& vAngle, Vec3& vPos)
{
	Vec3 vForward; Math::AngleVectors(vAngle, &vForward);
	const Vec3 vOrigin = pLocal->GetShootPos() + vForward * 128.f;

	CBaseEntity* pTarget;
	for (CEntitySphereQuery sphere(vOrigin, 128.f);
		(pTarget = sphere.GetCurrentEntity()) != nullptr;
		sphere.NextEntity())
	{
		if (pTarget == pEntity)
			break;
	}

	return pTarget == pEntity && SDK::VisPos(pLocal, pEntity, pLocal->GetShootPos(), vPos);
}

void CAutoAirblast::Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	if (!(Vars::Aimbot::Projectile::AutoAirblast.Value & Vars::Aimbot::Projectile::AutoAirblastEnum::Enabled) || !G::CanSecondaryAttack /*|| Vars::Auto::Airblast::DisableOnAttack.Value && pCmd->buttons & IN_ATTACK*/)
		return;

	const int iWeaponID = pWeapon->GetWeaponID();
	if (iWeaponID != TF_WEAPON_FLAMETHROWER && iWeaponID != TF_WEAPON_FLAME_BALL || SDK::AttribHookValue(0, "airblast_disabled", pWeapon))
		return;

	const Vec3 vEyePos = pLocal->GetShootPos();
	bool bShouldBlast = false;

	float flLatency = std::max(F::Backtrack.GetReal() - 0.05f, 0.f);

	for (auto pProjectile : H::Entities.GetGroup(EGroupType::WORLD_PROJECTILES))
	{
		if (!ShouldTargetProjectile(pProjectile, pLocal))
			continue;

		Vec3 vRestoreOrigin = pProjectile->GetAbsOrigin();
		Vec3 vOrigin = PredictOrigin(pProjectile->m_vecOrigin(), pProjectile->GetAbsVelocity(), flLatency);
		pProjectile->SetAbsOrigin(vOrigin);

		if (!(Vars::Aimbot::Projectile::AutoAirblast.Value & Vars::Aimbot::Projectile::AutoAirblastEnum::IgnoreFOV)
			&& Math::CalcFov(I::EngineClient->GetViewAngles(), Math::CalcAngle(vEyePos, vOrigin)) > Vars::Aimbot::General::AimFOV.Value)
			continue;

		if (Vars::Aimbot::Projectile::AutoAirblast.Value & Vars::Aimbot::Projectile::AutoAirblastEnum::Redirect) // implement advanced redirection
		{
			Vec3 vAngle = Math::CalcAngle(vEyePos, vOrigin);
			if (CanAirblastEntity(pLocal, pProjectile, vAngle, vOrigin))
			{
				SDK::FixMovement(pCmd, vAngle);
				pCmd->viewangles = vAngle;
				G::PSilentAngles = true;
				bShouldBlast = true;
			}
		}
		else if (CanAirblastEntity(pLocal, pProjectile, pCmd->viewangles, vOrigin))
			bShouldBlast = true;

		pProjectile->SetAbsOrigin(vRestoreOrigin);

		if (bShouldBlast)
			break;
	}

	if (bShouldBlast)
	{
		G::Attacking = true;
		pCmd->buttons |= IN_ATTACK2;
	}
}