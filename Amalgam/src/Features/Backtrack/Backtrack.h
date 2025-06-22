#pragma once
#include "../../SDK/SDK.h"

class CIncomingSequence
{
public:
	int m_nInReliableState;
	int m_nSequenceNr;
	float m_flTime;

	constexpr CIncomingSequence(int iState, int iSequence, float flTime) noexcept
		: m_nInReliableState{iState}, m_nSequenceNr{iSequence}, m_flTime{flTime}
	{
	}

	// Rule of 5: explicitly defaulted special member functions for better performance
	CIncomingSequence(const CIncomingSequence&) = default;
	CIncomingSequence& operator=(const CIncomingSequence&) = default;
	CIncomingSequence(CIncomingSequence&&) = default;
	CIncomingSequence& operator=(CIncomingSequence&&) = default;
	~CIncomingSequence() = default;
};

struct BoneMatrix
{
	matrix3x4 m_aBones[MAXSTUDIOBONES];
};

struct TickRecord
{
	float m_flSimTime = 0.0f;
	Vec3 m_vOrigin{};
	Vec3 m_vMins{};
	Vec3 m_vMaxs{};
	BoneMatrix m_BoneMatrix{};
	bool m_bOnShot = false;
	Vec3 m_vBreak{};
	bool m_bInvalid = false;

	// Constructor for better initialization
	TickRecord() noexcept = default;
	
	TickRecord(float simTime, const Vec3& origin, const Vec3& mins,
	                    const Vec3& maxs, const BoneMatrix& boneMatrix,
	                    bool onShot, const Vec3& breakVec) noexcept
		: m_flSimTime{simTime}, m_vOrigin{origin}, m_vMins{mins}, m_vMaxs{maxs},
		  m_BoneMatrix{boneMatrix}, m_bOnShot{onShot}, m_vBreak{breakVec}, m_bInvalid{false}
	{
	}

	// Rule of 5: explicitly defaulted for performance
	TickRecord(const TickRecord&) = default;
	TickRecord& operator=(const TickRecord&) = default;
	TickRecord(TickRecord&&) = default;
	TickRecord& operator=(TickRecord&&) = default;
	~TickRecord() = default;
};

class CBacktrack
{
private:
	// Private helper methods
	void UpdateDatagram();
	void MakeRecords();
	void CleanRecords();

	// Core data structures - using modern initialization
	std::unordered_map<CBaseEntity*, std::deque<TickRecord>> m_mRecords{};
	std::unordered_map<int, bool> m_mDidShoot{};

	// Sequence tracking
	std::deque<CIncomingSequence> m_dSequences{};
	int m_iLastInSequence{0};
	int m_nOldInSequenceNr{0};
	int m_nOldInReliableState{0};
	int m_nLastInSequenceNr{0};
	int m_nOldTickBase{0};
	
	// Configuration
	float m_flMaxUnlag{1.0f};

public:
	// Core functionality
	void Store();
	void SendLerp();
	void Draw(CTFPlayer* pLocal);
	void Reset() noexcept;

	// Record management
	[[nodiscard]] bool GetRecords(CBaseEntity* pEntity, std::vector<TickRecord*>& vReturn);
	[[nodiscard]] std::vector<TickRecord*> GetValidRecords(std::vector<TickRecord*>& vRecords,
	                                                       CTFPlayer* pLocal = nullptr,
	                                                       bool bDistance = false,
	                                                       float flTimeMod = 0.0f);

	// Latency and interpolation getters
	[[nodiscard]] constexpr float GetLerp() const noexcept;
	[[nodiscard]] constexpr float GetFake() const noexcept;
	[[nodiscard]] float GetReal(int iFlow = MAX_FLOWS, bool bNoFake = true) const noexcept;
	[[nodiscard]] constexpr float GetFakeInterp() const noexcept;
	[[nodiscard]] int GetAnticipatedChoke(int iMethod = Vars::Aimbot::General::AimType.Value) const noexcept;
	
	// Event handling
	void SetLerp(IGameEvent* pEvent);
	void ResolverUpdate(CBaseEntity* pEntity);
	void ReportShot(int iIndex);
	
	// Network manipulation
	void AdjustPing(CNetChannel* netChannel);
	void RestorePing(CNetChannel* netChannel) noexcept;

	// Public state - kept for compatibility
	int m_iTickCount{0};
	float m_flFakeLatency{0.0f};
	float m_flFakeInterp{0.015f};
	float m_flWishInterp{-1.0f};
	TickRecord m_tRecord{}; // for temporary use
};

ADD_FEATURE(CBacktrack, Backtrack);