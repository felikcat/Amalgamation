// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

MAKE_HOOK(IInput_GetUserCmd, U::Memory.GetVirtual(I::Input, 8), CUserCmd*,
	void* rcx, int sequence_number)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::IInput_GetUserCmd[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, sequence_number);
#endif

	return &I::Input->GetCommands()[sequence_number % MULTIPLAYER_BACKUP];
}