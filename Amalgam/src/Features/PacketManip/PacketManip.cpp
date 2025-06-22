// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "PacketManip.h"

#include "../Visuals/FakeAngle/FakeAngle.h"
#include "../Ticks/Ticks.h"

static inline bool AntiAimCheck(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd)
{
	return F::AntiAim.YawOn() && F::AntiAim.ShouldRun(pLocal, pWeapon, pCmd) && !F::Ticks.m_bRecharge
		&& I::ClientState->chokedcommands < F::AntiAim.AntiAimTicks();
}

void CPacketManip::Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd, bool* pSendPacket)
{
	F::FakeAngle.bDrawChams = Vars::Fakelag::Fakelag.Value || F::AntiAim.AntiAimOn();

	*pSendPacket = true;
	F::FakeLag.Run(pLocal, pWeapon, pCmd, pSendPacket);
	if (AntiAimCheck(pLocal, pWeapon, pCmd))
		*pSendPacket = false;
}