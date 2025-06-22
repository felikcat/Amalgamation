// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include <Windows.h>
#include "Core/Core.h"
#include "Utils/CrashLog/CrashLog.h"

DWORD WINAPI MainThread(LPVOID lpParam)
{
	U::Core.Load();
	U::Core.Loop();
	U::Core.Unload();

	CrashLog::Unload();
	FreeLibraryAndExitThread(static_cast<HMODULE>(lpParam), EXIT_SUCCESS);
}

// Do not try to make this "safe" as PVS Studio suggests.
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
	if (fdwReason == DLL_PROCESS_ATTACH)
	{
		CrashLog::Initialize();

		if (const auto hMainThread = CreateThread(nullptr, 0, MainThread, hinstDLL, 0, nullptr))
			CloseHandle(hMainThread);
	}

	return TRUE;
}