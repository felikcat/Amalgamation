// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

#include "../Features/Visuals/Materials/Materials.h"

MAKE_HOOK(IMatSystemSurface_OnScreenSizeChanged, U::Memory.GetVirtual(I::MatSystemSurface, 111), void,
	void* rcx, int nOldWidth, int nOldHeight)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::IMatSystemSurface_OnScreenSizeChanged[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, nOldWidth, nOldHeight);
#endif

	CALL_ORIGINAL(rcx, nOldWidth, nOldHeight);

	H::Fonts.Reload(Vars::Menu::Scale[DEFAULT_BIND]);
	F::Materials.ReloadMaterials();
}