// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "../SDK/SDK.h"

MAKE_HOOK(CStudioRender_SetAlphaModulation, U::Memory.GetVirtual(I::StudioRender, 28), void,
	void* rcx, float flAlpha)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::CStudioRender_SetAlphaModulation[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, flAlpha);
#endif

	if (Vars::Visuals::World::Modulations.Value & Vars::Visuals::World::ModulationsEnum::Prop && G::DrawingProps && !(Vars::Visuals::UI::CleanScreenshots.Value && I::EngineClient->IsTakingScreenshot()))
		return CALL_ORIGINAL(rcx, float(Vars::Colors::PropModulation.Value.a) / 255.f * flAlpha);

	CALL_ORIGINAL(rcx, flAlpha);
}