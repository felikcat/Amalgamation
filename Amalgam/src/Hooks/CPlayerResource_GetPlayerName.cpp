// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

#include "../Features/Players/PlayerUtils.h"

MAKE_SIGNATURE(CPlayerResource_GetPlayerName, "client.dll", "48 89 5C 24 ? 56 48 83 EC ? 48 63 F2", 0x0);

MAKE_HOOK(CPlayerResource_GetPlayerName, S::CPlayerResource_GetPlayerName(), const char*,
	void* rcx, int iIndex)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::CPlayerResource_GetPlayerName[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, iIndex);
#endif

	return F::PlayerUtils.GetPlayerName(iIndex, CALL_ORIGINAL(rcx, iIndex));
}