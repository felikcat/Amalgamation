#include "../SDK/SDK.h"

#include "../Features/Resolver/Resolver.h"
#include "../Features/EnginePrediction/EnginePrediction.h"

MAKE_HOOK(CBaseEntity_EstimateAbsVelocity, S::CBaseEntity_EstimateAbsVelocity(), void,
	void* rcx, Vector& vel)
{
#ifdef DEBUG_HOOKS
	if (!Vars::Hooks::CBaseEntity_EstimateAbsVelocity[DEFAULT_BIND])
		return CALL_ORIGINAL(rcx, vel);
#endif

	CALL_ORIGINAL(rcx, vel);

	auto pPlayer = reinterpret_cast<CTFPlayer*>(rcx);
	if (pPlayer->IsPlayer())
	{
		if (pPlayer->entindex() == I::EngineClient->GetLocalPlayer())
			vel = pPlayer->m_vecVelocity();
		else
		{
			// full override makes movement look a bit too choppy
			vel.z = pPlayer->m_vecVelocity().z;

			if (pPlayer->IsOnGround() && vel.Length2DSqr() < 2.f)
			{
				bool bMinwalk;
				if (F::Resolver.GetAngles(pPlayer, nullptr, nullptr, &bMinwalk) && bMinwalk)
				{
					// Use historical velocity data instead of artificial values
					if (auto pAvgVel = H::Entities.GetAvgVelocity(pPlayer->entindex()))
					{
						vel.x = pAvgVel->x;
						vel.y = pAvgVel->y;
						// Keep original Z for ground players
					}
					else
					{
						// Fallback to minimal artificial velocity only if no historical data
						vel = { 1, 1 };
					}
				}
			}
		}
	}
}