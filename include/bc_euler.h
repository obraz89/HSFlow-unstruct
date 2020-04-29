#pragma once

#include "bc_common.h"

#include "flux_euler.h"

enum struct t_BCKindEuler{
	Inflow=0,
	Outflow,
	Wall, 
	Sym
};

class t_BCEulerBase : public t_BCDataFace {
public:
	t_BCEulerBase() = delete;
	t_BCEulerBase(const std::string& fam_name) :t_BCDataFace(fam_name) {}
	virtual void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const= 0;
	virtual t_BCKindEuler getKind() const = 0;
};

// Supersonic Inflow
class t_BCDataInflow :public t_BCEulerBase {
public:

	t_BCDataInflow(const std::string& sect) : t_BCEulerBase(sect) {}
	// implement TPugin
	std::string get_description() const { return std::string("bc inflow"); };
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	// implement t_BCEulerBase
	void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindEuler getKindStat() { return t_BCKindEuler::Inflow; }
	t_BCKindEuler getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_inflow"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// supersonic outflow, i.e. extrapolation
class t_BCDataOutFlow :public t_BCEulerBase {

public:
	t_BCDataOutFlow(const std::string& sect) :t_BCEulerBase(sect) {}
	// implement TPugin
	std::string get_description() const { return std::string("bc outflow"); };
	// implement t_BCEulerBase
	void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindEuler getKindStat() { return t_BCKindEuler::Outflow; }
	t_BCKindEuler getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_outflow"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// Euler wall
class t_BCDataEulerWall :public t_BCEulerBase {

	double TWallDim;

public:
	t_BCDataEulerWall(const std::string& sect) :t_BCEulerBase(sect) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc euler wall"); };
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
	void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindEuler getKindStat() { return t_BCKindEuler::Wall; }
	t_BCKindEuler getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_euler_wall"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

// Symmetry bc
class t_BCDataSym :public t_BCEulerBase {
public:
	t_BCDataSym() = delete;
	t_BCDataSym(const std::string& sect) :t_BCEulerBase(sect) { }
	// implement TPugin
	std::string get_description() const { return std::string("bc sym"); };
	// implement t_BCEulerBase
	void yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const;
	static t_BCKindEuler getKindStat() { return t_BCKindEuler::Sym; }
	t_BCKindEuler getKind() const { return getKindStat(); }
	// implement t_BCDataFace
	static std::string getBCKindNameStat() { return "bc_sym"; }
	std::string getBCKindName() const { return getBCKindNameStat(); }
};

/**
 * List of bc sets for euler eqs
 // read bc_list section from config dynamically
 // so use iniFile capabilities directly instead of plugin params
 */

class t_BCListEuler : public t_BCList {
protected:
	std::vector<t_BCEulerBase*>  _pBCs;

public:
	// implement TPugin
	t_BCListEuler() {}
	std::string get_name() const { return "bc_list"; };
	std::string get_description() const { return std::string("bc list"); };
	// implement BCList
	std::string getSupportedBCsStr() const;

	t_FaceBCID getBCID(std::string fam_name) const;
	const t_BCEulerBase* getBC(int BCID) const;
	t_BCKindEuler getKind(int BCID) const { return getBC(BCID)->getKind(); }

	bool getBCKindByFamName(const std::string& name, t_BCKindEuler& kind) const;
	void addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data);
	bool has(std::string sectionname) const;
	void init(std::string& ini_data, const std::string& spec) {

		t_BCList::init(ini_data, spec);

	}

	virtual ~t_BCListEuler() {
		for (auto elem : _pBCs) delete elem;
	};

};

extern t_BCListEuler G_BCListEuler;

