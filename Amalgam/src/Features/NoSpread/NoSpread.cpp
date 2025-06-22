// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "NoSpread.h"

#include "NoSpreadProjectile/NoSpreadProjectile.h"
#include "NoSpreadHitscan/NoSpreadHitscan.h"

bool CNoSpread::ShouldRun(CTFPlayer* pLocal, CTFWeaponBase* pWeapon)
{
	if (!Vars::Aimbot::General::NoSpread.Value
		|| !pLocal || !pWeapon || !pLocal->CanAttack())
		return false;

	return true;
}

void CNoSpread::Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	if (!ShouldRun(pLocal, pWeapon))
		return;

	F::NoSpreadHitscan.Run(pLocal, pWeapon, pCmd);
	F::NoSpreadProjectile.Run(pLocal, pWeapon, pCmd);
}