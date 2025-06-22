// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

MAKE_SIGNATURE(CTFPlayer_AvoidPlayers, "client.dll", "48 89 54 24 ? 55 41 56", 0x0);

MAKE_HOOK(CTFPlayer_AvoidPlayers, S::CTFPlayer_AvoidPlayers(), void,
	void* rcx, CUserCmd* pCmd)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::CTFPlayer_AvoidPlayers[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, pCmd);
#endif

	if (Vars::Misc::Movement::NoPush.Value)
		return;

	CALL_ORIGINAL(rcx, pCmd);
}