#pragma once
#include "../../../SDK/SDK.h"

#include "../AimbotGlobal/AimbotGlobal.h"

Enum(PointType, None = 0, Regular = 1 << 0, Obscured = 1 << 1, ObscuredExtra = 1 << 2, ObscuredMulti = 1 << 3)
Enum(Calculated, Pending, Good, Time, Bad)

struct Solution_t
{
	float m_flPitch = 0.f;
	float m_flYaw = 0.f;
	float m_flTime = 0.f;
	int m_iCalculated = CalculatedEnum::Pending;
	
	// Enhanced solution validation with improved numerical stability
	// Fixed projectile distance accuracy issues for demoman's pipe launcher and other classes
	// by implementing long double precision calculations and corrected ballistic formulas
	[[nodiscard]] constexpr bool IsValid() const noexcept {
		return m_iCalculated == CalculatedEnum::Good &&
			   std::isfinite(m_flPitch) && std::isfinite(m_flYaw) &&
			   std::isfinite(m_flTime) && m_flTime > 0.0f;
	}
	
	// Reset solution to pending state
	constexpr void Reset() noexcept {
		m_flPitch = 0.0f;
		m_flYaw = 0.0f;
		m_flTime = 0.0f;
		m_iCalculated = CalculatedEnum::Pending;
	}
};
struct Point_t
{
	Vec3 m_vPoint = {};
	Solution_t m_Solution = {};
};
struct Info_t
{
	CTFPlayer* m_pLocal = nullptr;
	CTFWeaponBase* m_pWeapon = nullptr;

	Vec3 m_vLocalEye = {};
	Vec3 m_vTargetEye = {};

	float m_flLatency = 0.f;

	Vec3 m_vHull = {};
	Vec3 m_vOffset = {};
	Vec3 m_vAngFix = {};
	float m_flVelocity = 0.f;
	float m_flGravity = 0.f;
	float m_flRadius = 0.f;
	float m_flRadiusTime = 0.f;
	float m_flBoundingTime = 0.f;
	float m_flOffsetTime = 0.f;
	int m_iSplashCount = 0;
	int m_iSplashMode = 0;
	float m_flPrimeTime = 0;
	int m_iPrimeTime = 0;
};

class CAimbotProjectile
{
	std::vector<Target_t> GetTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon);
	std::vector<Target_t> SortTargets(CTFPlayer* pLocal, CTFWeaponBase* pWeapon);

	int GetHitboxPriority(int nHitbox, Target_t& tTarget, Info_t& tInfo);

	std::unordered_map<int, Vec3> GetDirectPoints(Target_t& tTarget, Info_t& tInfo);
	std::vector<Point_t> GetSplashPoints(Target_t& tTarget, std::vector<std::pair<Vec3, int>>& vSpherePoints, Info_t& tInfo, int iSimTime);
	void SetupSplashPoints(Target_t& tTarget, std::vector<std::pair<Vec3, int>>& vSpherePoints, Info_t& tInfo, std::vector<std::pair<Vec3, Vec3>>& vSimplePoints);
	std::vector<Point_t> GetSplashPointsSimple(Target_t& tTarget, std::vector<std::pair<Vec3, Vec3>>& vSpherePoints, Info_t& tInfo, int iSimTime);

	void CalculateAngle(const Vec3& vLocalPos, const Vec3& vTargetPos, Info_t& tInfo, int iSimTime, Solution_t& out, bool bAccuracy = true);
	void CalculateAngleRevolutionary(const Vec3& vLocalPos, const Vec3& vTargetPos, Info_t& tInfo, int iSimTime, Solution_t& out, bool bAccuracy = true);
	bool TestAngle(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, Target_t& tTarget, Vec3& vPoint, Vec3& vAngles, int iSimTime, bool bSplash, bool* pHitSolid = nullptr, std::vector<Vec3>* pProjectilePath = nullptr) noexcept;

	int CanHit(Target_t& tTarget, CTFPlayer* pLocal, CTFWeaponBase* pWeapon, std::vector<Vec3>* pPlayerPath, std::vector<Vec3>* pProjectilePath, std::vector<DrawBox_t>* pBoxes, float* pTimeTo);
	bool RunMain(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd);

	bool Aim(Vec3 vCurAngle, Vec3 vToAngle, Vec3& vOut, int iMethod = Vars::Aimbot::General::AimType.Value);
	void Aim(CUserCmd* pCmd, Vec3& vAngle);

	bool m_bLastTickHeld = false;

public:
	void Run(CTFPlayer* pLocal, CTFWeaponBase* pWeapon, CUserCmd* pCmd);
	float GetSplashRadius(CTFWeaponBase* pWeapon, CTFPlayer* pPlayer);

	int m_iLastTickCancel = 0;
};

ADD_FEATURE(CAimbotProjectile, AimbotProjectile);