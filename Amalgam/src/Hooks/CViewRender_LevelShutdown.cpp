// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

#include "../Features/Spectate/Spectate.h"

MAKE_HOOK(CViewRender_LevelShutdown, U::Memory.GetVirtual(I::ViewRender, 2), void,
	void* rcx)
{
	F::Spectate.m_iIntendedTarget = -1;

	CALL_ORIGINAL(rcx);
}