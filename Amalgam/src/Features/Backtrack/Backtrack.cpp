// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
#include "Backtrack.h"

#include "../PacketManip/FakeLag/FakeLag.h"
#include "../Ticks/Ticks.h"

void CBacktrack::Reset()
{
	m_mRecords.clear();
	m_dSequences.clear();
	m_iLastInSequence = 0;
}



// Enhanced adaptive interpolation for better accuracy
float CBacktrack::GetLerp()
{
	float flBaseInterp = Vars::Backtrack::Interp.Value / 1000.f;
	
	// Adaptive interpolation based on network conditions
	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (pNetChan)
	{
		const float flPacketLoss = pNetChan->GetAvgLoss(FLOW_INCOMING);
		const float flJitter = pNetChan->GetAvgLatency(FLOW_INCOMING) - pNetChan->GetLatency(FLOW_INCOMING);
		
		// Increase interpolation for poor network conditions
		if (flPacketLoss > 0.02f) // 2% packet loss threshold
			flBaseInterp *= (1.0f + flPacketLoss * 2.0f);
		
		if (fabsf(flJitter) > 0.01f) // 10ms jitter threshold
			flBaseInterp *= (1.0f + fabsf(flJitter) * 50.0f);
	}
	
	if (Vars::Misc::Game::AntiCheatCompatibility.Value)
		return std::clamp(flBaseInterp, G::Lerp, 0.1f);

	return std::clamp(flBaseInterp, G::Lerp, m_flMaxUnlag);
}

// Returns the wish backtrack latency
float CBacktrack::GetFake()
{
	return std::clamp(Vars::Backtrack::Latency.Value / 1000.f, 0.f, m_flMaxUnlag);
}

// Returns the current real latency
float CBacktrack::GetReal(int iFlow, bool bNoFake)
{
	auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan)
		return 0.f;

	if (iFlow != MAX_FLOWS)
		return pNetChan->GetLatency(iFlow) - (bNoFake && iFlow == FLOW_INCOMING ? m_flFakeLatency : 0.f);
	return pNetChan->GetLatency(FLOW_INCOMING) + pNetChan->GetLatency(FLOW_OUTGOING) - (bNoFake ? m_flFakeLatency : 0.f);
}

// Returns the current fake interp
float CBacktrack::GetFakeInterp()
{
	if (Vars::Misc::Game::AntiCheatCompatibility.Value)
		return std::min(m_flFakeInterp, 0.1f);

	return m_flFakeInterp;
}

// Enhanced anticipated choke calculation with server performance compensation
int CBacktrack::GetAnticipatedChoke(int iMethod)
{
	int iAnticipatedChoke = 0;
	if (F::Ticks.CanChoke() && G::PrimaryWeaponType != EWeaponType::HITSCAN && Vars::Aimbot::General::AimType.Value == Vars::Aimbot::General::AimTypeEnum::Silent)
		iAnticipatedChoke = 1;
	if (F::FakeLag.m_iGoal && !Vars::Fakelag::UnchokeOnAttack.Value && F::Ticks.m_iShiftedTicks == F::Ticks.m_iShiftedGoal && !F::Ticks.m_bDoubletap && !F::Ticks.m_bSpeedhack)
		iAnticipatedChoke = F::FakeLag.m_iGoal - I::ClientState->chokedcommands;
	
	// Dynamic server performance compensation
	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (pNetChan)
	{
		// Adjust for server tick variance
		const float flServerLatency = pNetChan->GetLatency(FLOW_OUTGOING);
		const float flTickVariance = flServerLatency - (flServerLatency > 0.05f ? 0.05f : 0.0f);
		
		// Add compensation for unstable servers
		if (flTickVariance > 0.02f) // 20ms variance threshold
		{
			const int iCompensation = static_cast<int>(flTickVariance / TICK_INTERVAL);
			iAnticipatedChoke += std::min(iCompensation, 3); // Cap at 3 ticks
		}
	}
	
	return iAnticipatedChoke;
}

void CBacktrack::SendLerp()
{
	static Timer tTimer = {};
	if (!tTimer.Run(0.1f))
		return;

	float flTarget = GetLerp();
	if (m_flWishInterp != flTarget)
	{
		m_flWishInterp = flTarget;

		auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
		if (pNetChan && I::EngineClient->IsConnected())
		{
			NET_SetConVar tConvar1 = { "cl_interp", std::to_string(m_flWishInterp).c_str() };
			pNetChan->SendNetMsg(tConvar1);

			NET_SetConVar tConvar2 = { "cl_interp_ratio", "1" };
			pNetChan->SendNetMsg(tConvar2);

			NET_SetConVar tConvar3 = { "cl_interpolate", "1" };
			pNetChan->SendNetMsg(tConvar3);
		}
	}
}

// Manages cl_interp client value
void CBacktrack::SetLerp(IGameEvent* pEvent)
{
	if (I::EngineClient->GetPlayerForUserID(pEvent->GetInt("userid")) == I::EngineClient->GetLocalPlayer())
		m_flFakeInterp = m_flWishInterp;
}

void CBacktrack::UpdateDatagram()
{
	const auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
	if (!pNetChan)
		return;

	if (const auto pLocal = H::Entities.GetLocal())
		m_nOldTickBase = pLocal->m_nTickBase();

	if (pNetChan->m_nInSequenceNr > m_iLastInSequence)
	{
		m_iLastInSequence = pNetChan->m_nInSequenceNr;
		m_dSequences.emplace_front(pNetChan->m_nInReliableState, pNetChan->m_nInSequenceNr, I::GlobalVars->realtime);
	}

	constexpr size_t maxSequences = 67;
	if (m_dSequences.size() > maxSequences)
		m_dSequences.pop_back();
}



bool CBacktrack::GetRecords(CBaseEntity* pEntity, std::vector<TickRecord*>& vReturn)
{
	const auto it = m_mRecords.find(pEntity);
	if (it == m_mRecords.end())
		return false;

	const auto& vRecords = it->second;
	vReturn.reserve(vReturn.size() + vRecords.size());
	
	for (const auto& tRecord : vRecords)
		vReturn.push_back(const_cast<TickRecord*>(&tRecord));
	
	return true;
}

std::vector<TickRecord*> CBacktrack::GetValidRecords(std::vector<TickRecord*>& vRecords, CTFPlayer* pLocal, bool bDistance, float flTimeMod)
{
	if (vRecords.empty())
		return {};

	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan)
		return {};

	std::vector<TickRecord*> vReturn;
	vReturn.reserve(vRecords.size()); // Pre-allocate to avoid reallocations
	
	const float flCorrect = std::clamp(GetReal(MAX_FLOWS, false) + ROUND_TO_TICKS(GetFakeInterp()), 0.f, m_flMaxUnlag);
	const int iServerTick = m_iTickCount + GetAnticipatedChoke() + Vars::Backtrack::Offset.Value + TIME_TO_TICKS(GetReal(FLOW_OUTGOING));
	const float flWindowThreshold = Vars::Backtrack::Window.Value / 1000.f;

	// Enhanced delta calculation with accuracy weighting
	const auto calculateWeightedDelta = [flCorrect, iServerTick, flTimeMod](const TickRecord* pRecord) -> std::pair<float, float> {
		const float flRecordTime = pRecord->m_flSimTime + flTimeMod;
		const float flTimeDiff = flCorrect - TICKS_TO_TIME(iServerTick - TIME_TO_TICKS(flRecordTime));
		const float flDelta = fabsf(flTimeDiff);
		
		// Calculate accuracy weight based on record quality
		float flWeight = 1.0f;
		
		// Prefer records that are not marked as invalid (higher quality)
		if (!pRecord->m_bInvalid)
			flWeight *= 1.5f;
		
		// Prefer records with shot events for better hit registration
		if (pRecord->m_bOnShot)
			flWeight *= 2.0f;
		
		// Prefer more recent records (within reason)
		const float flAge = flCorrect - pRecord->m_flSimTime;
		if (flAge < 0.1f) // Records within 100ms
			flWeight *= 1.3f;
		
		return {flDelta / flWeight, flWeight}; // Return weighted delta and weight
	};

	if (!Vars::Misc::Game::AntiCheatCompatibility.Value && Vars::Backtrack::Window.Value)
	{
		// Enhanced filtering with quality consideration
		for (const auto pRecord : vRecords)
		{
			const auto [weightedDelta, weight] = calculateWeightedDelta(pRecord);
			if (weightedDelta <= flWindowThreshold)
				vReturn.push_back(pRecord);
		}
	}

	if (vReturn.empty())
	{	// Enhanced fallback selection with quality metrics
		constexpr float flMaxDelta = 0.2f;
		float flBestScore = std::numeric_limits<float>::max();
		TickRecord* pBestRecord = nullptr;
		
		for (const auto pRecord : vRecords)
		{
			const auto [weightedDelta, weight] = calculateWeightedDelta(pRecord);
			
			// Combined score considering both time accuracy and record quality
			const float flScore = weightedDelta + (pRecord->m_bInvalid ? 0.1f : 0.0f);
			
			if (flScore < flBestScore && weightedDelta <= flMaxDelta)
			{
				flBestScore = flScore;
				pBestRecord = pRecord;
			}
		}
		
		if (pBestRecord)
			vReturn.push_back(pBestRecord);
	}
	else if (pLocal && vReturn.size() > 1)
	{
		const auto localOrigin = pLocal->m_vecOrigin();
		const bool bPreferOnShot = Vars::Backtrack::PreferOnShot.Value;
		
		if (bDistance)
		{
			std::sort(vReturn.begin(), vReturn.end(), [bPreferOnShot, &localOrigin](const TickRecord* a, const TickRecord* b) -> bool
				{
					if (bPreferOnShot && a->m_bOnShot != b->m_bOnShot)
						return a->m_bOnShot > b->m_bOnShot;

					// Use squared distance to avoid expensive sqrt operations
					const float distSqrA = localOrigin.DistToSqr(a->m_vOrigin);
					const float distSqrB = localOrigin.DistToSqr(b->m_vOrigin);
					return distSqrA < distSqrB;
				});
		}
		else
		{
			std::sort(vReturn.begin(), vReturn.end(), [bPreferOnShot, calculateWeightedDelta](const TickRecord* a, const TickRecord* b) -> bool
				{
					if (bPreferOnShot && a->m_bOnShot != b->m_bOnShot)
						return a->m_bOnShot > b->m_bOnShot;

					const auto [deltaA, weightA] = calculateWeightedDelta(a);
					const auto [deltaB, weightB] = calculateWeightedDelta(b);
					return deltaA < deltaB;
				});
		}
	}

	return vReturn;
}



void CBacktrack::MakeRecords()
{
	static auto sv_lagcompensation_teleport_dist = U::ConVars.FindVar("sv_lagcompensation_teleport_dist");
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	for (const auto& pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		const int playerIndex = pPlayer->entindex();
		
		if (playerIndex == localPlayerIndex || pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost()
			|| !H::Entities.GetDeltaTime(playerIndex))
			continue;

		const auto pBones = H::Entities.GetBones(playerIndex);
		if (!pBones)
			continue;

		auto& vRecords = m_mRecords[pPlayer];
		const TickRecord* pLastRecord = vRecords.empty() ? nullptr : &vRecords.front();
		
		// Create new record efficiently
		vRecords.emplace_front(
			pPlayer->m_flSimulationTime(),
			pPlayer->m_vecOrigin(),
			pPlayer->m_vecMins(),
			pPlayer->m_vecMaxs(),
			*reinterpret_cast<const BoneMatrix*>(pBones),
			m_mDidShoot[playerIndex],
			pPlayer->m_vecOrigin()
		);
		const auto& tCurRecord = vRecords.front();

		// Enhanced accuracy features for backtracking
		if (pLastRecord)
		{
			// Velocity-based prediction for better accuracy
			const Vec3 vVelocity = (tCurRecord.m_vOrigin - pLastRecord->m_vOrigin) / (tCurRecord.m_flSimTime - pLastRecord->m_flSimTime);
			const float flSpeed = vVelocity.Length();
			
			// Detect suspicious movement patterns that might indicate cheating or network issues
			constexpr float maxReasonableSpeed = 1000.f; // Units per second
			const bool bSuspiciousMovement = flSpeed > maxReasonableSpeed;
			
			// Enhanced record validation based on movement consistency
			if (bSuspiciousMovement)
			{
				// Mark recent records as potentially inaccurate but don't invalidate completely
				for (auto& tRecord : vRecords)
				{
					if (tRecord.m_flSimTime > tCurRecord.m_flSimTime - 0.5f) // Last 500ms
						tRecord.m_bInvalid = true;
				}
			}
			
			// Adaptive record quality scoring
			float flRecordQuality = 1.0f;
			
			// Reduce quality for high-speed movement
			if (flSpeed > 400.f) // Typical running speed threshold
				flRecordQuality *= std::max(0.3f, 1.0f - (flSpeed - 400.f) / 600.f);
			
			// Improve quality for consistent movement patterns
			if (vRecords.size() >= 3)
			{
				const auto& secondLast = vRecords[1];
				const Vec3 vPrevVelocity = (pLastRecord->m_vOrigin - secondLast.m_vOrigin) / (pLastRecord->m_flSimTime - secondLast.m_flSimTime);
				const float flVelocityConsistency = 1.0f - std::min(1.0f, vVelocity.DistTo(vPrevVelocity) / 500.f);
				flRecordQuality *= (0.7f + 0.3f * flVelocityConsistency);
			}
			
			// Store quality metric in the current record (using unused field creatively)
			const_cast<TickRecord&>(tCurRecord).m_bInvalid = flRecordQuality < 0.5f;
			
			// Predictive backtracking enhancement for improved accuracy
			if (vRecords.size() >= 2)
			{
				// Calculate movement prediction for next tick
				const float flTimeDelta = tCurRecord.m_flSimTime - pLastRecord->m_flSimTime;
				if (flTimeDelta > 0.0f && flTimeDelta < 0.1f) // Reasonable time delta
				{
					// Store predicted position for improved accuracy with fast-moving targets
					const float flPredictionTime = TICK_INTERVAL;
					const Vec3 vPredictedPos = tCurRecord.m_vOrigin + vVelocity * flPredictionTime;
					
					// Store prediction in break position for efficiency
					const_cast<TickRecord&>(tCurRecord).m_vBreak = vPredictedPos;
				}
			}
		}
		m_mDidShoot[playerIndex] = false;
	}
}

void CBacktrack::CleanRecords()
{
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	// Pre-calculate deadtime threshold to avoid repeated calculations
	const float flCurrentTime = I::GlobalVars->curtime;
	const float flRealLatency = GetReal();
	const float flDeadtime = flCurrentTime + flRealLatency - m_flMaxUnlag;
	
	constexpr float flMaxSimTime = std::numeric_limits<float>::max();
	
	for (const auto& pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		if (pPlayer->entindex() == localPlayerIndex)
			continue;

		auto& vRecords = m_mRecords[pPlayer];

		if (pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost())
		{
			vRecords.clear();
			continue;
		}

		// Remove hack record if it exists and there are other records
		if (vRecords.size() > 1 && vRecords.back().m_flSimTime == flMaxSimTime)
			vRecords.pop_back();
		
		// Remove expired records from the back
		while (!vRecords.empty())
		{
			const auto& backRecord = vRecords.back();
			if (backRecord.m_flSimTime < flDeadtime ||
				(vRecords.size() > 1 && backRecord.m_flSimTime == flMaxSimTime))
			{
				vRecords.pop_back();
			}
			else
			{
				break;
			}
		}
	}
}



void CBacktrack::Store()
{
	UpdateDatagram();
	if (!I::EngineClient->IsInGame())
		return;

	static const auto sv_maxunlag = U::ConVars.FindVar("sv_maxunlag");
	m_flMaxUnlag = sv_maxunlag->GetFloat();
	
	MakeRecords();
	CleanRecords();
}

void CBacktrack::ResolverUpdate(CBaseEntity* pEntity)
{
	/*
	if (!m_mRecords.contains(pEntity))
		return;

	m_mRecords[pEntity].clear();
	*/
}

void CBacktrack::ReportShot(int iIndex)
{
	if (!Vars::Backtrack::PreferOnShot.Value)
		return;

	const auto pEntity = I::ClientEntityList->GetClientEntity(iIndex);
	if (!pEntity)
		return;

	const auto pPlayer = pEntity->As<CTFPlayer>();
	const auto pWeapon = pPlayer->m_hActiveWeapon().Get();
	if (!pWeapon || SDK::GetWeaponType(pWeapon->As<CTFWeaponBase>()) != EWeaponType::HITSCAN)
		return;

	m_mDidShoot[pEntity->entindex()] = true;
}

void CBacktrack::AdjustPing(CNetChannel* pNetChan)
{
	m_nOldInSequenceNr = pNetChan->m_nInSequenceNr;
	m_nOldInReliableState = pNetChan->m_nInReliableState;

	const auto calculateLatency = [this, pNetChan]() -> float
		{
			if (!Vars::Backtrack::Latency.Value)
				return 0.f;

			const auto pLocal = H::Entities.GetLocal();
			if (!pLocal || !pLocal->m_iClass())
				return 0.f;

			static const auto host_timescale = U::ConVars.FindVar("host_timescale");
			const float flTimescale = host_timescale->GetFloat();

			static float flStaticReal = 0.f;
			const float flFake = GetFake();
			const float flReal = TICKS_TO_TIME(pLocal->m_nTickBase() - m_nOldTickBase);
			
			// Use compile-time constants for better optimization
			constexpr float smoothingFactor = 0.1f;
			constexpr float tickMultiplier = 5.0f;
			const float tickOffset = tickMultiplier * TICK_INTERVAL;
			
			// Optimize the smoothing calculation
			const float flTarget = flReal + tickOffset;
			flStaticReal += (flTarget - flStaticReal) * smoothingFactor;

			int nInReliableState = pNetChan->m_nInReliableState;
			int nInSequenceNr = pNetChan->m_nInSequenceNr;
			float flLatency = 0.f;
			
			const float flMaxLatency = m_flMaxUnlag - flStaticReal;
			
			// Cache realtime to avoid repeated global access
			const float flRealtime = I::GlobalVars->realtime;
			
			for (const auto& cSequence : m_dSequences)
			{
				nInReliableState = cSequence.m_nInReliableState;
				nInSequenceNr = cSequence.m_nSequenceNr;
				
				// Optimize latency calculation by reducing operations
				const float flTimeDiff = flRealtime - cSequence.m_flTime;
				flLatency = flTimeDiff * flTimescale - TICK_INTERVAL;

				if (flLatency > flFake || m_nLastInSequenceNr >= cSequence.m_nSequenceNr || flLatency > flMaxLatency)
					break;
			}

			pNetChan->m_nInReliableState = nInReliableState;
			pNetChan->m_nInSequenceNr = nInSequenceNr;
			return flLatency;
		};

	const float flLatency = calculateLatency();
	m_nLastInSequenceNr = pNetChan->m_nInSequenceNr;

	if (Vars::Backtrack::Latency.Value || m_flFakeLatency)
	{
		constexpr float smoothingRate = 0.1f;
		const float flDelta = (flLatency - m_flFakeLatency) * smoothingRate;
		const float flNewValue = m_flFakeLatency + flDelta;
		
		// Pre-calculate clamp bounds to avoid repeated arithmetic
		const float flLowerBound = m_flFakeLatency - TICK_INTERVAL;
		const float flUpperBound = m_flFakeLatency + TICK_INTERVAL;
		
		m_flFakeLatency = std::clamp(flNewValue, flLowerBound, flUpperBound);
		
		// Use epsilon comparison for better floating-point handling
		constexpr float epsilon = 1e-6f;
		if (flLatency < epsilon && m_flFakeLatency < TICK_INTERVAL)
			m_flFakeLatency = 0.f;
	}
}

void CBacktrack::RestorePing(CNetChannel* pNetChan)
{
	pNetChan->m_nInSequenceNr = m_nOldInSequenceNr;
	pNetChan->m_nInReliableState = m_nOldInReliableState;
}

void CBacktrack::Draw(CTFPlayer* pLocal)
{
	if (!(Vars::Menu::Indicators.Value & Vars::Menu::IndicatorsEnum::Ping) || !pLocal->IsAlive())
		return;

	const auto pResource = H::Entities.GetPR();
	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pResource || !pNetChan)
		return;

	static float flFakeLatency = 0.f;
	{
		static Timer tTimer = {};
		if (tTimer.Run(0.5f))
			flFakeLatency = m_flFakeLatency;
	}
	
	const float flFakeInterp = GetFakeInterp();
	const float flFakeLerp = flFakeInterp > G::Lerp ? flFakeInterp : 0.f;

	// Use compile-time constant for millisecond conversion
	constexpr float msConversion = 1000.f;
	
	// Optimize calculations by reducing function calls and intermediate operations
	const float flFakeTotal = flFakeLatency + flFakeLerp;
	const float flFake = std::min(flFakeTotal, m_flMaxUnlag) * msConversion;
	
	const float flTotalLatency = pNetChan->GetLatency(FLOW_INCOMING) + pNetChan->GetLatency(FLOW_OUTGOING);
	const float flLatency = std::max(flTotalLatency - flFakeLatency, 0.f) * msConversion;
	const int iLatencyScoreboard = pResource->m_iPing(pLocal->entindex());

	int x = Vars::Menu::PingDisplay.Value.x;
	int y = Vars::Menu::PingDisplay.Value.y + 8;
	const auto& fFont = H::Fonts.GetFont(FONT_INDICATORS);
	const int nTall = fFont.m_nTall + H::Draw.Scale(1);

	// Determine alignment based on screen position
	EAlign align = ALIGN_TOP;
	constexpr int edgeThreshold = 100;
	const int scaledOffset = H::Draw.Scale(50, Scale_Round);
	const int adjustmentOffset = H::Draw.Scale(42, Scale_Round);
	
	if (x <= edgeThreshold + scaledOffset)
	{
		x -= adjustmentOffset;
		align = ALIGN_TOPLEFT;
	}
	else if (x >= H::Draw.m_nScreenW - edgeThreshold - scaledOffset)
	{
		x += adjustmentOffset;
		align = ALIGN_TOPRIGHT;
	}

	// Draw ping information
	constexpr size_t bufferSize = 64;
	char szBuffer[bufferSize];
	
	const auto& activeColor = Vars::Menu::Theme::Active.Value;
	const auto& bgColor = Vars::Menu::Theme::Background.Value;
	
	if (flFake > 0.f || Vars::Backtrack::Interp.Value)
	{
		sprintf_s(szBuffer, bufferSize, "Ping %.0f (+ %.0f) ms", flLatency, flFake);
	}
	else
	{
		sprintf_s(szBuffer, bufferSize, "Ping %.0f ms", flLatency);
	}
	H::Draw.StringOutlined(fFont, x, y, activeColor, bgColor, align, szBuffer);
	
	// Draw scoreboard ping
	y += nTall;
	sprintf_s(szBuffer, bufferSize, "Scoreboard %d ms", iLatencyScoreboard);
	H::Draw.StringOutlined(fFont, x, y, activeColor, bgColor, align, szBuffer);
}