#include "Backtrack.h"

#include "../PacketManip/FakeLag/FakeLag.h"
#include "../Ticks/Ticks.h"

// Ensure std::format is available
#include <format>

void CBacktrack::Reset()
{
	m_mRecords.clear();
	m_dSequences.clear();
	m_iLastInSequence = 0;
	m_nOldInSequenceNr = 0;
	m_nOldInReliableState = 0;
	m_nLastInSequenceNr = 0;
	m_nOldTickBase = 0;
	m_flFakeLatency = 0.0f;
	m_flFakeInterp = 0.015f;
	m_flWishInterp = -1.0f;
	m_iTickCount = 0;
	m_tRecord = {};
}



// Returns the wish cl_interp
float CBacktrack::GetLerp()
{
	const float interpValue = Vars::Backtrack::Interp.Value / 1000.0f;
	
	if (Vars::Misc::Game::AntiCheatCompatibility.Value)
		return std::clamp(interpValue, G::Lerp, 0.1f);

	return std::clamp(interpValue, G::Lerp, m_flMaxUnlag);
}

// Returns the wish backtrack latency
float CBacktrack::GetFake()
{
	const float latencyValue = Vars::Backtrack::Latency.Value / 1000.0f;
	return std::clamp(latencyValue, 0.0f, m_flMaxUnlag);
}

// Returns the current real latency
float CBacktrack::GetReal(int iFlow, bool bNoFake)
{
	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan)
		return 0.0f;

	if (iFlow != MAX_FLOWS)
	{
		const float latency = pNetChan->GetLatency(iFlow);
		const float fakeAdjustment = (bNoFake && iFlow == FLOW_INCOMING) ? m_flFakeLatency : 0.0f;
		return latency - fakeAdjustment;
	}
	
	const float incomingLatency = pNetChan->GetLatency(FLOW_INCOMING);
	const float outgoingLatency = pNetChan->GetLatency(FLOW_OUTGOING);
	const float fakeAdjustment = bNoFake ? m_flFakeLatency : 0.0f;
	return incomingLatency + outgoingLatency - fakeAdjustment;
}

// Returns the current fake interp
float CBacktrack::GetFakeInterp()
{
	constexpr float maxCompatibilityInterp = 0.1f;
	
	if (Vars::Misc::Game::AntiCheatCompatibility.Value)
		return std::min(m_flFakeInterp, maxCompatibilityInterp);

	return m_flFakeInterp;
}

// Returns the current anticipated choke
int CBacktrack::GetAnticipatedChoke(int iMethod)
{
	int anticipatedChoke = 0;
	
	// Check for silent aim choke
	const bool canSilentChoke = F::Ticks.CanChoke() &&
	                           G::PrimaryWeaponType != EWeaponType::HITSCAN &&
	                           Vars::Aimbot::General::AimType.Value == Vars::Aimbot::General::AimTypeEnum::Silent;
	if (canSilentChoke)
		anticipatedChoke = 1;
	
	// Check for fake lag choke
	const bool canFakeLagChoke = F::FakeLag.m_iGoal &&
	                            !Vars::Fakelag::UnchokeOnAttack.Value &&
	                            F::Ticks.m_iShiftedTicks == F::Ticks.m_iShiftedGoal &&
	                            !F::Ticks.m_bDoubletap &&
	                            !F::Ticks.m_bSpeedhack;
	if (canFakeLagChoke)
		anticipatedChoke = F::FakeLag.m_iGoal - I::ClientState->chokedcommands;
	
	return anticipatedChoke;
}

void CBacktrack::SendLerp()
{
	static Timer timer{};
	constexpr float updateInterval = 0.1f;
	
	if (!timer.Run(updateInterval))
		return;

	const float targetInterp = GetLerp();
	if (m_flWishInterp != targetInterp)
	{
		m_flWishInterp = targetInterp;

		const auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
		if (pNetChan && I::EngineClient->IsConnected())
		{
			// Send interpolation settings
			NET_SetConVar interpConVar{"cl_interp", std::to_string(m_flWishInterp).c_str()};
			pNetChan->SendNetMsg(interpConVar);

			NET_SetConVar ratioConVar{"cl_interp_ratio", "1"};
			pNetChan->SendNetMsg(ratioConVar);

			NET_SetConVar enableConVar{"cl_interpolate", "1"};
			pNetChan->SendNetMsg(enableConVar);
		}
	}
}

// Manages cl_interp client value
void CBacktrack::SetLerp(IGameEvent* pEvent)
{
	if (!pEvent)
		return;
		
	const int userID = pEvent->GetInt("userid");
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	if (I::EngineClient->GetPlayerForUserID(userID) == localPlayerIndex)
		m_flFakeInterp = m_flWishInterp;
}

void CBacktrack::UpdateDatagram()
{
	const auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
	if (!pNetChan)
		return;

	// Update tick base from local player
	if (const auto pLocal = H::Entities.GetLocal())
		m_nOldTickBase = pLocal->m_nTickBase();

	// Add new sequence if it's newer
	if (pNetChan->m_nInSequenceNr > m_iLastInSequence)
	{
		m_iLastInSequence = pNetChan->m_nInSequenceNr;
		m_dSequences.emplace_front(
			pNetChan->m_nInReliableState,
			pNetChan->m_nInSequenceNr,
			I::GlobalVars->realtime
		);
	}

	// Limit sequence history size
	constexpr size_t maxSequences = 67;
	if (m_dSequences.size() > maxSequences)
		m_dSequences.pop_back();
}



bool CBacktrack::GetRecords(CBaseEntity* pEntity, std::vector<TickRecord*>& vReturn)
{
	const auto it = m_mRecords.find(pEntity);
	if (it == m_mRecords.end())
		return false;

	auto& vRecords = it->second;
	vReturn.reserve(vRecords.size()); // Pre-allocate for better performance
	for (auto& tRecord : vRecords)
		vReturn.push_back(&tRecord);
	return true;
}

std::vector<TickRecord*> CBacktrack::GetValidRecords(std::vector<TickRecord*>& vRecords, CTFPlayer* pLocal, bool bDistance, float flTimeMod)
{
	if (vRecords.empty())
		return {};

	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan)
		return {};

	std::vector<TickRecord*> vReturn{};
	vReturn.reserve(vRecords.size()); // Pre-allocate for better performance
	
	const float flCorrect = std::clamp(GetReal(MAX_FLOWS, false) + ROUND_TO_TICKS(GetFakeInterp()), 0.0f, m_flMaxUnlag);
	const int iServerTick = m_iTickCount + GetAnticipatedChoke() + Vars::Backtrack::Offset.Value + TIME_TO_TICKS(GetReal(FLOW_OUTGOING));
	const float flWindowThreshold = Vars::Backtrack::Window.Value / 1000.0f;

	// Lambda to calculate delta for reuse
	const auto calculateDelta = [&](const TickRecord* pRecord) -> float {
		return std::abs(flCorrect - TICKS_TO_TIME(iServerTick - TIME_TO_TICKS(pRecord->m_flSimTime + flTimeMod)));
	};

	if (!Vars::Misc::Game::AntiCheatCompatibility.Value && Vars::Backtrack::Window.Value)
	{
		for (const auto* pRecord : vRecords)
		{
			const float flDelta = calculateDelta(pRecord);
			if (flDelta <= flWindowThreshold)
				vReturn.push_back(const_cast<TickRecord*>(pRecord));
		}
	}

	if (vReturn.empty())
	{	// Ensure there is at least 1 record
		constexpr float flMaxDelta = 0.2f;
		float flMinDelta = flMaxDelta;
		TickRecord* pBestRecord = nullptr;
		
		for (auto* pRecord : vRecords)
		{
			const float flDelta = calculateDelta(pRecord);
			if (flDelta <= flMinDelta)
			{
				flMinDelta = flDelta;
				pBestRecord = pRecord;
			}
		}
		
		if (pBestRecord)
			vReturn.push_back(pBestRecord);
	}
	else if (pLocal && vReturn.size() > 1)
	{
		// Sort records based on preference
		if (bDistance)
		{
			std::sort(vReturn.begin(), vReturn.end(), [&](const TickRecord* a, const TickRecord* b) -> bool
			{
				if (Vars::Backtrack::PreferOnShot.Value && a->m_bOnShot != b->m_bOnShot)
					return a->m_bOnShot > b->m_bOnShot;

				return pLocal->m_vecOrigin().DistTo(a->m_vOrigin) < pLocal->m_vecOrigin().DistTo(b->m_vOrigin);
			});
		}
		else
		{
			std::sort(vReturn.begin(), vReturn.end(), [&](const TickRecord* a, const TickRecord* b) -> bool
			{
				if (Vars::Backtrack::PreferOnShot.Value && a->m_bOnShot != b->m_bOnShot)
					return a->m_bOnShot > b->m_bOnShot;

				const float flADelta = calculateDelta(a);
				const float flBDelta = calculateDelta(b);
				return flADelta < flBDelta;
			});
		}
	}

	return vReturn;
}



void CBacktrack::MakeRecords()
{
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	for (const auto& pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		const int playerIndex = pPlayer->entindex();
		
		// Early exit conditions
		if (playerIndex == localPlayerIndex || pPlayer->IsDormant() ||
		    !pPlayer->IsAlive() || pPlayer->IsAGhost() ||
		    !H::Entities.GetDeltaTime(playerIndex))
			continue;

		const auto pBones = H::Entities.GetBones(playerIndex);
		if (!pBones)
			continue;

		auto& vRecords = m_mRecords[pPlayer];
		const TickRecord* pLastRecord = vRecords.empty() ? nullptr : &vRecords.front();
		
		// Create new record using constructor
		vRecords.emplace_front(
			pPlayer->m_flSimulationTime(),
			pPlayer->m_vecOrigin(),
			pPlayer->m_vecMins(),
			pPlayer->m_vecMaxs(),
			*reinterpret_cast<const BoneMatrix*>(pBones),
			m_mDidShoot[playerIndex],
			pPlayer->m_vecOrigin()
		);
		
		const TickRecord& currentRecord = vRecords.front();
		bool bLagComp = false;
		
		if (pLastRecord)
		{
			const Vec3 deltaPosition = currentRecord.m_vBreak - pLastRecord->m_vBreak;
			
			static const auto sv_lagcompensation_teleport_dist = U::ConVars.FindVar("sv_lagcompensation_teleport_dist");
			const float teleportDistSqr = std::pow(sv_lagcompensation_teleport_dist->GetFloat(), 2.0f);
			
			if (deltaPosition.Length2DSqr() > teleportDistSqr)
			{
				bLagComp = true;
				if (!H::Entities.GetLagCompensation(playerIndex))
				{
					vRecords.resize(1);
					vRecords.front().m_flSimTime = std::numeric_limits<float>::max(); // Invalidation marker
				}
				
				// Mark all records as invalid using range-based for loop
				for (auto& record : vRecords)
					record.m_bInvalid = true;
			}

			// Update invalid records with current data
			for (auto& record : vRecords)
			{
				if (record.m_bInvalid)
				{
					record.m_vOrigin = currentRecord.m_vOrigin;
					record.m_vMins = currentRecord.m_vMins;
					record.m_vMaxs = currentRecord.m_vMaxs;
					record.m_BoneMatrix = currentRecord.m_BoneMatrix;
					record.m_bOnShot = currentRecord.m_bOnShot;
				}
			}
		}

		H::Entities.SetLagCompensation(playerIndex, bLagComp);
		m_mDidShoot[playerIndex] = false;
	}
}

void CBacktrack::CleanRecords()
{
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	const float deadTime = I::GlobalVars->curtime + GetReal() - m_flMaxUnlag;
	constexpr float invalidTime = std::numeric_limits<float>::max();

	for (const auto& pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		if (pPlayer->entindex() == localPlayerIndex)
			continue;

		auto& vRecords = m_mRecords[pPlayer];

		// Clear records for invalid players
		if (pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost())
		{
			vRecords.clear();
			continue;
		}

		// Remove invalid marker record if it exists
		if (vRecords.size() > 1 && vRecords.back().m_flSimTime == invalidTime)
			vRecords.pop_back();

		// Remove expired records
		while (!vRecords.empty())
		{
			const auto& backRecord = vRecords.back();
			const bool isExpired = backRecord.m_flSimTime < deadTime;
			const bool isInvalidMarker = (vRecords.size() > 1) && (backRecord.m_flSimTime == invalidTime);
			
			if (isExpired || isInvalidMarker)
				vRecords.pop_back();
			else
				break;
		}
	}
}



void CBacktrack::Store()
{
	UpdateDatagram();
	if (!I::EngineClient->IsInGame())
		return;

	static const auto sv_maxunlag = U::ConVars.FindVar("sv_maxunlag");
	if (sv_maxunlag)
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
	if (!pPlayer)
		return;
		
	const auto pWeapon = pPlayer->m_hActiveWeapon().Get();
	if (!pWeapon)
		return;
		
	if (SDK::GetWeaponType(pWeapon->As<CTFWeaponBase>()) == EWeaponType::HITSCAN)
		m_mDidShoot[pEntity->entindex()] = true;
}

void CBacktrack::AdjustPing(CNetChannel* pNetChan)
{
	// Store original values
	m_nOldInSequenceNr = pNetChan->m_nInSequenceNr;
	m_nOldInReliableState = pNetChan->m_nInReliableState;

	const auto calculateLatency = [this, pNetChan]() -> float
	{
		if (!Vars::Backtrack::Latency.Value)
			return 0.0f;

		const auto pLocal = H::Entities.GetLocal();
		if (!pLocal || !pLocal->m_iClass())
			return 0.0f;

		static const auto host_timescale = U::ConVars.FindVar("host_timescale");
		const float timescale = host_timescale ? host_timescale->GetFloat() : 1.0f;

		static float staticReal = 0.0f;
		const float fakeLatency = GetFake();
		const float realLatency = TICKS_TO_TIME(pLocal->m_nTickBase() - m_nOldTickBase);
		
		// Smooth the static real value
		constexpr float smoothingFactor = 0.1f;
		const float baseOffset = 5.0f * TICK_INTERVAL;
		staticReal += (realLatency + baseOffset - staticReal) * smoothingFactor;

		int bestReliableState = pNetChan->m_nInReliableState;
		int bestSequenceNr = pNetChan->m_nInSequenceNr;
		float bestLatency = 0.0f;

		for (const auto& sequence : m_dSequences)
		{
			const float currentLatency = (I::GlobalVars->realtime - sequence.m_flTime) * timescale - TICK_INTERVAL;
			
			// Check if this sequence is valid
			if (currentLatency > fakeLatency ||
			    m_nLastInSequenceNr >= sequence.m_nSequenceNr ||
			    currentLatency > m_flMaxUnlag - staticReal)
				break;

			bestReliableState = sequence.m_nInReliableState;
			bestSequenceNr = sequence.m_nSequenceNr;
			bestLatency = currentLatency;
		}

		pNetChan->m_nInReliableState = bestReliableState;
		pNetChan->m_nInSequenceNr = bestSequenceNr;
		return bestLatency;
	};

	const float calculatedLatency = calculateLatency();
	m_nLastInSequenceNr = pNetChan->m_nInSequenceNr;

	// Update fake latency with smoothing
	if (Vars::Backtrack::Latency.Value || m_flFakeLatency > 0.0f)
	{
		constexpr float smoothingRate = 0.1f;
		const float deltaLatency = (calculatedLatency - m_flFakeLatency) * smoothingRate;
		const float clampedDelta = std::clamp(deltaLatency, -TICK_INTERVAL, TICK_INTERVAL);
		
		m_flFakeLatency = std::clamp(m_flFakeLatency + clampedDelta, 0.0f, m_flMaxUnlag);
		
		// Reset if latency is very small
		if (calculatedLatency == 0.0f && m_flFakeLatency < TICK_INTERVAL)
			m_flFakeLatency = 0.0f;
	}
}

void CBacktrack::RestorePing(CNetChannel* pNetChan)
{
	if (!pNetChan)
		return;
		
	pNetChan->m_nInSequenceNr = m_nOldInSequenceNr;
	pNetChan->m_nInReliableState = m_nOldInReliableState;
}

void CBacktrack::Draw(CTFPlayer* pLocal)
{
	if (!(Vars::Menu::Indicators.Value & Vars::Menu::IndicatorsEnum::Ping) || !pLocal->IsAlive())
		return;

	auto pResource = H::Entities.GetPR();
	auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pResource || !pNetChan)
		return;

	static float flFakeLatency = 0.f;
	{
		static Timer tTimer = {};
		if (tTimer.Run(0.5f))
			flFakeLatency = m_flFakeLatency;
	}
	float flFakeLerp = GetFakeInterp() > G::Lerp ? GetFakeInterp() : 0.f;

	float flFake = std::min(flFakeLatency + flFakeLerp, m_flMaxUnlag) * 1000.f;
	float flLatency = std::max(pNetChan->GetLatency(FLOW_INCOMING) + pNetChan->GetLatency(FLOW_OUTGOING) - flFakeLatency, 0.f) * 1000.f;
	int iLatencyScoreboard = pResource->m_iPing(pLocal->entindex());

	int x = Vars::Menu::PingDisplay.Value.x;
	int y = Vars::Menu::PingDisplay.Value.y + 8;
	const auto& fFont = H::Fonts.GetFont(FONT_INDICATORS);
	const int nTall = fFont.m_nTall + H::Draw.Scale(1);

	EAlign align = ALIGN_TOP;
	if (x <= 100 + H::Draw.Scale(50, Scale_Round))
	{
		x -= H::Draw.Scale(42, Scale_Round);
		align = ALIGN_TOPLEFT;
	}
	else if (x >= H::Draw.m_nScreenW - 100 - H::Draw.Scale(50, Scale_Round))
	{
		x += H::Draw.Scale(42, Scale_Round);
		align = ALIGN_TOPRIGHT;
	}

	if (flFake || Vars::Backtrack::Interp.Value)
		H::Draw.StringOutlined(fFont, x, y, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, std::format("Ping {:.0f} (+ {:.0f}) ms", flLatency, flFake).c_str());
	else
		H::Draw.StringOutlined(fFont, x, y, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, std::format("Ping {:.0f} ms", flLatency).c_str());
	H::Draw.StringOutlined(fFont, x, y += nTall, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, std::format("Scoreboard {} ms", iLatencyScoreboard).c_str());
}