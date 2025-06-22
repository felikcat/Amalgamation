#pragma once
#include "../../../SDK/SDK.h"

Enum(ESPText, Top, Bottom, Right, Health, Uber)

struct PlayerCache
{
	float m_flAlpha = 1.f;
	std::vector<std::tuple<int, std::string, Color_t, Color_t>> m_vText;
	Color_t m_tColor = {};
	bool m_bBox = false;
	bool m_bBones = false;

	bool m_bHealthBar = false;
	bool m_bUberBar = false;
	int m_iClassIcon = 0;
	CHudTexture* m_pWeaponIcon = nullptr;
	float m_flHealth = 1.f;
	float m_flUber = 0.f;

	// Reserve space for common text entries to avoid reallocations
	void ReserveText() noexcept { m_vText.reserve(8); }
	void ClearText() noexcept { m_vText.clear(); }
};

struct BuildingCache
{
	float m_flAlpha = 1.f;
	std::vector<std::tuple<int, std::string, Color_t, Color_t>> m_vText;
	Color_t m_tColor = {};
	bool m_bBox = false;

	bool m_bHealthBar = false;
	float m_flHealth = 1.f;

	// Reserve space for common text entries to avoid reallocations
	void ReserveText() noexcept { m_vText.reserve(4); }
	void ClearText() noexcept { m_vText.clear(); }
};

struct WorldCache
{
	float m_flAlpha = 1.f;
	std::vector<std::tuple<int, std::string, Color_t, Color_t>> m_vText;
	Color_t m_tColor = {};
	bool m_bBox = false;

	// Reserve space for common text entries to avoid reallocations
	void ReserveText() noexcept { m_vText.reserve(2); }
	void ClearText() noexcept { m_vText.clear(); }
};

class CESP
{
private:
	void StorePlayers(CTFPlayer* pLocal);
	void StoreBuildings(CTFPlayer* pLocal);
	void StoreProjectiles(CTFPlayer* pLocal);
	void StoreObjective(CTFPlayer* pLocal);
	void StoreWorld();

	void DrawPlayers();
	void DrawBuildings();
	void DrawWorld();

	Color_t GetColor(CTFPlayer* pLocal, CBaseEntity* pEntity);
	bool GetDrawBounds(CBaseEntity* pEntity, float& x, float& y, float& w, float& h);
	void DrawBones(CTFPlayer* pPlayer, matrix3x4* aBones, std::vector<int> vecBones, Color_t clr);

	// Inline helper for better performance
	inline const char* GetPlayerClass(int nClassNum) const noexcept
	{
		static constexpr const char* szClasses[] = {
			"Unknown", "Scout", "Sniper", "Soldier", "Demoman", "Medic", "Heavy", "Pyro", "Spy", "Engineer"
		};
		return (nClassNum < 10 && nClassNum > 0) ? szClasses[nClassNum] : szClasses[0];
	}

	// Reserve initial capacity to reduce allocations
	std::unordered_map<CBaseEntity*, PlayerCache> m_mPlayerCache;
	std::unordered_map<CBaseEntity*, BuildingCache> m_mBuildingCache;
	std::unordered_map<CBaseEntity*, WorldCache> m_mWorldCache;

public:
	CESP() noexcept
	{
		// Reserve space for common entity counts to reduce hash map reallocations
		m_mPlayerCache.reserve(32);    // Typical max players in TF2
		m_mBuildingCache.reserve(16);  // Reasonable building count
		m_mWorldCache.reserve(64);     // World entities can be numerous
	}

	void Store(CTFPlayer* pLocal);
	void Draw();
};

ADD_FEATURE(CESP, ESP);