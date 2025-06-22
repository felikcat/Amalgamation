// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include <Windows.h>
#include <process.h>
#include "Core/Core.h"
#include "Utils/CrashLog/CrashLog.h"

unsigned __stdcall MainThread(void* lpParam)
{
	U::Core.Load();
	U::Core.Loop();
	U::Core.Unload();

	CrashLog::Unload();
	FreeLibraryAndExitThread(static_cast<HMODULE>(lpParam), EXIT_SUCCESS);
	return 0;
}

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
	if (fdwReason == DLL_PROCESS_ATTACH)
	{
		CrashLog::Initialize();

		// Use _beginthreadex instead of CreateThread to avoid the PVS-Studio V718 warning:
		// _beginthreadex is safer for use with CRT functions and doesn't hold the loader lock.
		if (const auto hMainThread = reinterpret_cast<HANDLE>(_beginthreadex(nullptr, 0, MainThread, hinstDLL, 0, nullptr)))
			CloseHandle(hMainThread);
	}

	return TRUE;
}