// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "Timer.h"

#include "../../SDK/SDK.h"

Timer::Timer()
{
	m_flLast = SDK::PlatFloatTime();
}

bool Timer::Check(float flS) const
{
	float flCurrentTime = SDK::PlatFloatTime();
	return flCurrentTime - m_flLast >= flS;
}

bool Timer::Run(float flS)
{
	if (Check(flS))
	{
		Update();
		return true;
	}
	return false;
}

inline void Timer::Update()
{
	m_flLast = SDK::PlatFloatTime();
}