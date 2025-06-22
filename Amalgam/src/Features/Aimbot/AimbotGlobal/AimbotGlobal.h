#pragma once
#include "../../../SDK/SDK.h"
#include "../../Backtrack/Backtrack.h"

Enum(Target, Unknown, Player, Sentry, Dispenser, Teleporter, Sticky, NPC, Bomb)

struct Target_t
{
	CBaseEntity* m_pEntity = nullptr;
	int m_iTargetType = TargetEnum::Unknown;
	Vec3 m_vPos = {};
	Vec3 m_vAngleTo = {};
	float m_flFOVTo = std::numeric_limits<float>::max();
	float m_flDistTo = std::numeric_limits<float>::max();
	int m_nPriority = 0;
	int m_nAimedHitbox = -1;

	TickRecord* m_pRecord = nullptr;
	bool m_bBacktrack = false;
};

class CAimbotGlobal
{
public:
	void SortTargets(std::vector<Target_t>&, const int iMethod) noexcept;
	void SortPriority(std::vector<Target_t>&) noexcept;

	bool PlayerBoneInFOV(CTFPlayer* pTarget, const Vec3 vLocalPos, const Vec3 vLocalAngles, float& flFOVTo, Vec3& vPos, Vec3& vAngleTo, const int iHitboxes = Vars::Aimbot::Hitscan::HitboxesEnum::Head | Vars::Aimbot::Hitscan::HitboxesEnum::Body | Vars::Aimbot::Hitscan::HitboxesEnum::Pelvis | Vars::Aimbot::Hitscan::HitboxesEnum::Arms | Vars::Aimbot::Hitscan::HitboxesEnum::Legs) noexcept;
	bool IsHitboxValid(const uint32_t uHash, const int nHitbox, const int iHitboxes = Vars::Aimbot::Hitscan::HitboxesEnum::Head | Vars::Aimbot::Hitscan::HitboxesEnum::Body | Vars::Aimbot::Hitscan::HitboxesEnum::Pelvis | Vars::Aimbot::Hitscan::HitboxesEnum::Arms | Vars::Aimbot::Hitscan::HitboxesEnum::Legs) noexcept;

	bool ShouldIgnore(CBaseEntity* pTarget, CTFPlayer* pLocal, CTFWeaponBase* pWeapon) noexcept;
	int GetPriority(const int iIndex) noexcept;

	bool ShouldAim() noexcept;
	bool ShouldHoldAttack(CTFWeaponBase* pWeapon) noexcept;
	bool ValidBomb(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CBaseEntity* pBomb) noexcept;
};

ADD_FEATURE(CAimbotGlobal, AimbotGlobal);