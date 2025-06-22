// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

MAKE_SIGNATURE(CInventoryManager_ShowItemsPickedUp, "client.dll", "44 88 4C 24 ? 44 88 44 24 ? 53 41 56", 0x0);

MAKE_HOOK(CInventoryManager_ShowItemsPickedUp, S::CInventoryManager_ShowItemsPickedUp(), bool,
	void* rcx, bool bForce, bool bReturnToGame, bool bNoPanel)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::CInventoryManager_ShowItemsPickedUp[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, bForce, bReturnToGame, bNoPanel);
#endif

	if (Vars::Misc::Automation::AcceptItemDrops.Value)
	{
		CALL_ORIGINAL(rcx, true, true, true);
		return false;
	}
	return CALL_ORIGINAL(rcx, bForce, bReturnToGame, bNoPanel);
}