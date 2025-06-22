#include "Backtrack.h"

#include "../PacketManip/FakeLag/FakeLag.h"
#include "../Ticks/Ticks.h"

#include <algorithm>
#include <numeric>
#include <ranges>
#include <concepts>
#include <span>
#include <expected>
#include <sstream>
#include <iomanip>

void CBacktrack::Reset() noexcept
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
[[nodiscard]] constexpr float CBacktrack::GetLerp() const noexcept
{
	constexpr float millisToSeconds = 0.001f;
	constexpr float maxCompatibilityInterp = 0.1f;
	
	const float interpValue = Vars::Backtrack::Interp.Value * millisToSeconds;
	
	return Vars::Misc::Game::AntiCheatCompatibility.Value
		? std::clamp(interpValue, G::Lerp, maxCompatibilityInterp)
		: std::clamp(interpValue, G::Lerp, m_flMaxUnlag);
}

// Returns the wish backtrack latency
[[nodiscard]] constexpr float CBacktrack::GetFake() const noexcept
{
	constexpr float millisToSeconds = 0.001f;
	const float latencyValue = Vars::Backtrack::Latency.Value * millisToSeconds;
	return std::clamp(latencyValue, 0.0f, m_flMaxUnlag);
}

// Returns the current real latency
[[nodiscard]] float CBacktrack::GetReal(int iFlow, bool bNoFake) const noexcept
{
	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan) [[unlikely]]
		return 0.0f;

	if (iFlow != MAX_FLOWS) [[likely]]
	{
		const float latency = pNetChan->GetLatency(iFlow);
		const float fakeAdjustment = (bNoFake && iFlow == FLOW_INCOMING) ? m_flFakeLatency : 0.0f;
		return std::max(0.0f, latency - fakeAdjustment);
	}
	
	const float incomingLatency = pNetChan->GetLatency(FLOW_INCOMING);
	const float outgoingLatency = pNetChan->GetLatency(FLOW_OUTGOING);
	const float fakeAdjustment = bNoFake ? m_flFakeLatency : 0.0f;
	return std::max(0.0f, incomingLatency + outgoingLatency - fakeAdjustment);
}

// Returns the current fake interp
[[nodiscard]] constexpr float CBacktrack::GetFakeInterp() const noexcept
{
	constexpr float maxCompatibilityInterp = 0.1f;
	
	return Vars::Misc::Game::AntiCheatCompatibility.Value
		? std::min(m_flFakeInterp, maxCompatibilityInterp)
		: m_flFakeInterp;
}

// Returns the current anticipated choke
[[nodiscard]] int CBacktrack::GetAnticipatedChoke(int iMethod) const noexcept
{
	// Check for silent aim choke
	const bool canSilentChoke = F::Ticks.CanChoke() &&
	                           G::PrimaryWeaponType != EWeaponType::HITSCAN &&
	                           Vars::Aimbot::General::AimType.Value == Vars::Aimbot::General::AimTypeEnum::Silent;
	if (canSilentChoke) [[unlikely]]
		return 1;
	
	// Check for fake lag choke
	const bool canFakeLagChoke = F::FakeLag.m_iGoal &&
	                            !Vars::Fakelag::UnchokeOnAttack.Value &&
	                            F::Ticks.m_iShiftedTicks == F::Ticks.m_iShiftedGoal &&
	                            !F::Ticks.m_bDoubletap &&
	                            !F::Ticks.m_bSpeedhack;
	if (canFakeLagChoke) [[likely]]
		return std::max(0, F::FakeLag.m_iGoal - I::ClientState->chokedcommands);
	
	return 0;
}

void CBacktrack::SendLerp()
{
	static Timer timer{};
	constexpr float updateInterval = 0.1f;
	
	if (!timer.Run(updateInterval)) [[likely]]
		return;

	const float targetInterp = GetLerp();
	if (m_flWishInterp == targetInterp) [[likely]]
		return;
		
	m_flWishInterp = targetInterp;

	const auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
	if (!pNetChan || !I::EngineClient->IsConnected()) [[unlikely]]
		return;

	// Use structured bindings and modern initialization for better readability
	const auto sendConVar = [pNetChan](std::string_view name, std::string_view value) {
		NET_SetConVar conVar{name.data(), value.data()};
		pNetChan->SendNetMsg(conVar);
	};

	// Send interpolation settings using the lambda
	sendConVar("cl_interp", std::to_string(m_flWishInterp));
	sendConVar("cl_interp_ratio", "1");
	sendConVar("cl_interpolate", "1");
}

// Manages cl_interp client value
void CBacktrack::SetLerp(IGameEvent* pEvent)
{
	if (!pEvent) [[unlikely]]
		return;
		
	const int userID = pEvent->GetInt("userid");
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	if (I::EngineClient->GetPlayerForUserID(userID) == localPlayerIndex) [[unlikely]]
		m_flFakeInterp = m_flWishInterp;
}

void CBacktrack::UpdateDatagram()
{
	const auto pNetChan = reinterpret_cast<CNetChannel*>(I::EngineClient->GetNetChannelInfo());
	if (!pNetChan) [[unlikely]]
		return;

	// Update tick base from local player using structured binding
	if (const auto pLocal = H::Entities.GetLocal(); pLocal) [[likely]]
		m_nOldTickBase = pLocal->m_nTickBase();

	// Add new sequence if it's newer
	if (pNetChan->m_nInSequenceNr > m_iLastInSequence) [[likely]]
	{
		m_iLastInSequence = pNetChan->m_nInSequenceNr;
		m_dSequences.emplace_front(
			pNetChan->m_nInReliableState,
			pNetChan->m_nInSequenceNr,
			I::GlobalVars->realtime
		);

		// Limit sequence history size immediately after adding
		constexpr size_t maxSequences = 67;
		if (m_dSequences.size() > maxSequences) [[unlikely]]
			m_dSequences.pop_back();
	}
}



bool CBacktrack::GetRecords(CBaseEntity* pEntity, std::vector<TickRecord*>& vReturn)
{
	if (const auto it = m_mRecords.find(pEntity); it != m_mRecords.end()) [[likely]]
	{
		auto& vRecords = it->second;
		vReturn.reserve(vRecords.size()); // Pre-allocate for better performance
		
		// Use transform for more expressive code
		std::transform(vRecords.begin(), vRecords.end(), std::back_inserter(vReturn),
			[](auto& record) -> TickRecord* { return &record; });
		return true;
	}
	return false;
}

std::vector<TickRecord*> CBacktrack::GetValidRecords(std::vector<TickRecord*>& vRecords, CTFPlayer* pLocal, bool bDistance, float flTimeMod)
{
	if (vRecords.empty()) [[unlikely]]
		return {};

	const auto pNetChan = I::EngineClient->GetNetChannelInfo();
	if (!pNetChan) [[unlikely]]
		return {};

	std::vector<TickRecord*> vReturn{};
	vReturn.reserve(vRecords.size()); // Pre-allocate for better performance
	
	constexpr float millisToSeconds = 0.001f;
	const float flCorrect = std::clamp(GetReal(MAX_FLOWS, false) + ROUND_TO_TICKS(GetFakeInterp()), 0.0f, m_flMaxUnlag);
	const int iServerTick = m_iTickCount + GetAnticipatedChoke() + Vars::Backtrack::Offset.Value + TIME_TO_TICKS(GetReal(FLOW_OUTGOING));
	const float flWindowThreshold = Vars::Backtrack::Window.Value * millisToSeconds;

	// Lambda to calculate delta for reuse
	const auto calculateDelta = [&](const TickRecord* pRecord) -> float {
		return std::abs(flCorrect - TICKS_TO_TIME(iServerTick - TIME_TO_TICKS(pRecord->m_flSimTime + flTimeMod)));
	};

	// Use copy_if for more expressive filtering
	if (!Vars::Misc::Game::AntiCheatCompatibility.Value && Vars::Backtrack::Window.Value) [[likely]]
	{
		std::copy_if(vRecords.begin(), vRecords.end(), std::back_inserter(vReturn),
			[&](const TickRecord* pRecord) {
				return calculateDelta(pRecord) <= flWindowThreshold;
			});
	}

	if (vReturn.empty()) [[unlikely]]
	{
		// Find best record using min_element
		constexpr float flMaxDelta = 0.2f;
		
		if (const auto bestIt = std::min_element(vRecords.begin(), vRecords.end(),
			[&](const TickRecord* a, const TickRecord* b) {
				return calculateDelta(a) < calculateDelta(b);
			}); bestIt != vRecords.end() && calculateDelta(*bestIt) <= flMaxDelta)
		{
			vReturn.push_back(*bestIt);
		}
	}
	else if (pLocal && vReturn.size() > 1) [[likely]]
	{
		// Sort records based on preference using modern lambda with capture
		const auto createComparator = [&](bool useDistance) {
			return [&, useDistance](const TickRecord* a, const TickRecord* b) -> bool {
				if (Vars::Backtrack::PreferOnShot.Value && a->m_bOnShot != b->m_bOnShot)
					return a->m_bOnShot > b->m_bOnShot;

				return useDistance
					? pLocal->m_vecOrigin().DistTo(a->m_vOrigin) < pLocal->m_vecOrigin().DistTo(b->m_vOrigin)
					: calculateDelta(a) < calculateDelta(b);
			};
		};

		std::sort(vReturn.begin(), vReturn.end(), createComparator(bDistance));
	}

	return vReturn;
}



void CBacktrack::MakeRecords()
{
	const int localPlayerIndex = I::EngineClient->GetLocalPlayer();
	
	// Use structured binding and modern range-based for loop
	for (const auto& pEntity : H::Entities.GetGroup(EGroupType::PLAYERS_ALL))
	{
		const auto pPlayer = pEntity->As<CTFPlayer>();
		const int playerIndex = pPlayer->entindex();
		
		// Early exit conditions with [[likely]] hints
		if (playerIndex == localPlayerIndex || pPlayer->IsDormant() ||
		    !pPlayer->IsAlive() || pPlayer->IsAGhost() ||
		    !H::Entities.GetDeltaTime(playerIndex)) [[likely]]
			continue;

		const auto pBones = H::Entities.GetBones(playerIndex);
		if (!pBones) [[unlikely]]
			continue;

		auto& vRecords = m_mRecords[pPlayer];
		const TickRecord* pLastRecord = vRecords.empty() ? nullptr : &vRecords.front();
		
		// Create new record using constructor with perfect forwarding
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
		
		if (pLastRecord) [[likely]]
		{
			// Teleport detection removed - allow all movement patterns
			bLagComp = false;

			// Update invalid records with current data using algorithm
			const auto updateInvalidRecord = [&currentRecord](auto& record) {
				if (record.m_bInvalid) [[unlikely]]
				{
					record.m_vOrigin = currentRecord.m_vOrigin;
					record.m_vMins = currentRecord.m_vMins;
					record.m_vMaxs = currentRecord.m_vMaxs;
					record.m_BoneMatrix = currentRecord.m_BoneMatrix;
					record.m_bOnShot = currentRecord.m_bOnShot;
				}
			};
			
			std::for_each(vRecords.begin(), vRecords.end(), updateInvalidRecord);
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
		if (pPlayer->entindex() == localPlayerIndex) [[unlikely]]
			continue;

		auto& vRecords = m_mRecords[pPlayer];

		// Clear records for invalid players
		if (pPlayer->IsDormant() || !pPlayer->IsAlive() || pPlayer->IsAGhost()) [[unlikely]]
		{
			vRecords.clear();
			continue;
		}

		// Remove invalid marker record if it exists
		if (vRecords.size() > 1 && vRecords.back().m_flSimTime == invalidTime) [[unlikely]]
			vRecords.pop_back();

		// Remove expired records using modern predicate-based approach
		const auto isExpiredOrInvalid = [&](const TickRecord& record) -> bool {
			const bool isExpired = record.m_flSimTime < deadTime;
			const bool isInvalidMarker = (vRecords.size() > 1) && (record.m_flSimTime == invalidTime);
			return isExpired || isInvalidMarker;
		};

		// Use reverse iterator to remove from back efficiently
		auto reverseIt = std::find_if_not(vRecords.rbegin(), vRecords.rend(), isExpiredOrInvalid);
		if (reverseIt != vRecords.rbegin()) [[likely]]
		{
			vRecords.erase(reverseIt.base(), vRecords.end());
		}
	}
}



void CBacktrack::Store()
{
	UpdateDatagram();
	if (!I::EngineClient->IsInGame()) [[unlikely]]
		return;

	// Cache the ConVar lookup for better performance
	static const auto sv_maxunlag = U::ConVars.FindVar("sv_maxunlag");
	if (sv_maxunlag) [[likely]]
		m_flMaxUnlag = sv_maxunlag->GetFloat();
	
	MakeRecords();
	CleanRecords();
}

void CBacktrack::ResolverUpdate(CBaseEntity* pEntity)
{
	// Check if entity exists in records and clear if found
	if (const auto it = m_mRecords.find(pEntity); it != m_mRecords.end()) [[unlikely]]
	{
		it->second.clear();
	}
}

void CBacktrack::ReportShot(int iIndex)
{
	if (!Vars::Backtrack::PreferOnShot.Value) [[likely]]
		return;

	const auto pEntity = I::ClientEntityList->GetClientEntity(iIndex);
	if (!pEntity) [[unlikely]]
		return;
		
	const auto pPlayer = pEntity->As<CTFPlayer>();
	if (!pPlayer) [[unlikely]]
		return;
		
	const auto pWeapon = pPlayer->m_hActiveWeapon().Get();
	if (!pWeapon) [[unlikely]]
		return;
		
	if (SDK::GetWeaponType(pWeapon->As<CTFWeaponBase>()) == EWeaponType::HITSCAN) [[likely]]
		m_mDidShoot[pEntity->entindex()] = true;
}

void CBacktrack::AdjustPing(CNetChannel* pNetChan)
{
	// Store original values using structured binding
	const auto [oldSequence, oldReliable] = std::make_pair(pNetChan->m_nInSequenceNr, pNetChan->m_nInReliableState);
	m_nOldInSequenceNr = oldSequence;
	m_nOldInReliableState = oldReliable;

	const auto calculateLatency = [this, pNetChan]() -> float
	{
		if (!Vars::Backtrack::Latency.Value) [[likely]]
			return 0.0f;

		const auto pLocal = H::Entities.GetLocal();
		if (!pLocal || !pLocal->m_iClass()) [[unlikely]]
			return 0.0f;

		static const auto host_timescale = U::ConVars.FindVar("host_timescale");
		const float timescale = host_timescale ? host_timescale->GetFloat() : 1.0f;

		static float staticReal = 0.0f;
		const float fakeLatency = GetFake();
		const float realLatency = TICKS_TO_TIME(pLocal->m_nTickBase() - m_nOldTickBase);
		
		// Smooth the static real value using modern constexpr
		constexpr float smoothingFactor = 0.1f;
		const float baseOffset = 5.0f * TICK_INTERVAL;
		staticReal += (realLatency + baseOffset - staticReal) * smoothingFactor;

		// Use structured binding for better readability
		auto [bestReliableState, bestSequenceNr, bestLatency] =
			std::make_tuple(pNetChan->m_nInReliableState, pNetChan->m_nInSequenceNr, 0.0f);

		// Use modern algorithm with early termination
		const auto isValidSequence = [&](const auto& sequence) -> bool {
			const float currentLatency = (I::GlobalVars->realtime - sequence.m_flTime) * timescale - TICK_INTERVAL;
			
			return !(currentLatency > fakeLatency ||
			         m_nLastInSequenceNr >= sequence.m_nSequenceNr ||
			         currentLatency > m_flMaxUnlag - staticReal);
		};

		for (const auto& sequence : m_dSequences)
		{
			const float currentLatency = (I::GlobalVars->realtime - sequence.m_flTime) * timescale - TICK_INTERVAL;
			
			if (!isValidSequence(sequence)) [[unlikely]]
				break;

			std::tie(bestReliableState, bestSequenceNr, bestLatency) =
				std::make_tuple(sequence.m_nInReliableState, sequence.m_nSequenceNr, currentLatency);
		}

		std::tie(pNetChan->m_nInReliableState, pNetChan->m_nInSequenceNr) =
			std::make_pair(bestReliableState, bestSequenceNr);
		return bestLatency;
	};

	const float calculatedLatency = calculateLatency();
	m_nLastInSequenceNr = pNetChan->m_nInSequenceNr;

	// Update fake latency with smoothing using modern approach
	if (Vars::Backtrack::Latency.Value || m_flFakeLatency > 0.0f) [[likely]]
	{
		constexpr float smoothingRate = 0.1f;
		const float deltaLatency = (calculatedLatency - m_flFakeLatency) * smoothingRate;
		const float clampedDelta = std::clamp(deltaLatency, -TICK_INTERVAL, TICK_INTERVAL);
		
		m_flFakeLatency = std::clamp(m_flFakeLatency + clampedDelta, 0.0f, m_flMaxUnlag);
		
		// Reset if latency is very small using epsilon comparison
		constexpr float epsilon = 1e-6f;
		if (std::abs(calculatedLatency) < epsilon && m_flFakeLatency < TICK_INTERVAL) [[unlikely]]
			m_flFakeLatency = 0.0f;
	}
}

void CBacktrack::RestorePing(CNetChannel* pNetChan) noexcept
{
	if (!pNetChan) [[unlikely]]
		return;
		
	// Use structured binding for cleaner restoration
	std::tie(pNetChan->m_nInSequenceNr, pNetChan->m_nInReliableState) =
		std::make_pair(m_nOldInSequenceNr, m_nOldInReliableState);
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

	// Helper lambda for string formatting
	const auto formatString = [](const auto&... args) -> std::string {
		std::ostringstream oss;
		oss << std::fixed << std::setprecision(0);
		((oss << args), ...);
		return oss.str();
	};

	if (flFake || Vars::Backtrack::Interp.Value) {
		const auto pingText = formatString("Ping ", flLatency, " (+ ", flFake, ") ms");
		H::Draw.StringOutlined(fFont, x, y, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, pingText.c_str());
	} else {
		const auto pingText = formatString("Ping ", flLatency, " ms");
		H::Draw.StringOutlined(fFont, x, y, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, pingText.c_str());
	}
	
	const auto scoreboardText = formatString("Scoreboard ", iLatencyScoreboard, " ms");
	H::Draw.StringOutlined(fFont, x, y += nTall, Vars::Menu::Theme::Active.Value, Vars::Menu::Theme::Background.Value, align, scoreboardText.c_str());
}