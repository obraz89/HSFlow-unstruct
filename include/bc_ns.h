#pragma once

#pragma once

#include "bc_common.h"

#include "flux_ns.h"

enum t_BCKindNS {
	InflowSup = 0,
	OutflowSup,
	WallNoSlip,
	EulerWall
};

class t_BCNSBase : public t_BCDataFace {
public:
	t_BCNSBase() = delete;
	t_BCNSBase(const std::string& fam_name) :t_BCDataFace(fam_name) {}
	//virtual void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const = 0;
	virtual t_BCKindNS getKind() const = 0;
};

// Supersonic Inflow
class t_BCNSInflowSup :public t_BCNSBase {
public:

	t_BCNSInflowSup(const std::string& sect) : t_BCNSBase(sect) {}
	// implement TPugin
	std::string get_description() const { return std::string("bc inflow"); };
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	// implement t_BCEulerBase
	//void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindNS getKindStat() { return t_BCKindNS::InflowSup; }
	t_BCKindNS getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_ns_inflow_sup"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// supersonic outflow, i.e. extrapolation
class t_BCNSOutflowSup :public t_BCNSBase {

public:
	t_BCNSOutflowSup(const std::string& sect) :t_BCNSBase(sect) {}
	// implement TPugin
	std::string get_description() const { return std::string("bc outflow"); };
	// implement t_BCEulerBase
	//void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindNS getKindStat() { return t_BCKindNS::OutflowSup; }
	t_BCKindNS getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_ns_outflow_sup"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// No-slip wall
class t_BCNSWall :public t_BCNSBase {

	double TWallDim;

public:
	t_BCNSWall(const std::string& sect) :t_BCNSBase(sect) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc ns noslip wall"); };
	void default_settings() {
		TPluginParamsGroup g("", "gas-dynamic functions values on the wall");
		g.add("Tw_K", 300.0, "temperature, dimensional in K");
		_mapParamsGrps.emplace(g.get_name(), g);
	};
	void init(std::string& ini_data, const std::string& spec) {

		TPlugin::init(ini_data, spec);

		const TPluginParamsGroup& g = get_settings_grp("");

		TWallDim = g.get_real_param("Tw_K");

	};
	// implement t_BCEulerBase
	//void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindNS getKindStat() { return t_BCKindNS::WallNoSlip; }
	t_BCKindNS getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_ns_wall"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// Euler Wall
class t_BCNSEuWall :public t_BCNSBase {
public:
	t_BCNSEuWall() = delete;
	t_BCNSEuWall(const std::string& sect) :t_BCNSBase(sect) { }
	// implement TPugin
	std::string get_description() const { return std::string("bc ns euler wall"); };
	// implement t_BCNSBase
	//void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindNS getKindStat() { return t_BCKindNS::EulerWall; }
	t_BCKindNS getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_ns_euler_wall"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

/**
 * List of bc sets for euler eqs
 // read bc_list section from config dynamically
 // so use iniFile capabilities directly instead of plugin params
 */

class t_BCListNS : public t_BCList {
protected:
	std::vector<t_BCNSBase*>  _pBCs;

public:
	// implement TPugin
	t_BCListNS() {}
	std::string get_name() const { return "bc_list"; };
	std::string get_description() const { return std::string("bc list"); };
	// implement BCList
	std::string getSupportedBCsStr() const;

	t_FaceBCID getBCID(std::string fam_name) const;
	const t_BCNSBase* getBC(int BCID) const;
	t_BCKindNS getKind(int BCID) const { return getBC(BCID)->getKind(); }

	bool getBCKindByFamName(const std::string& name, t_BCKindNS& kind) const;
	void addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data);
	bool has(std::string sectionname) const;
	void init(std::string& ini_data, const std::string& spec) {

		t_BCList::init(ini_data, spec);

	}

	virtual ~t_BCListNS() {
		for (auto elem : _pBCs) delete elem;
	};

};

extern t_BCListNS G_BCListNS;


