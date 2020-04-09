#pragma once

#include "bc_common.h"

#include "flux_euler.h"

// This is BCId, always start from 1 as 0 is for fluid faces
enum struct t_BCKindEuler{
	Inflow=0,
	Outflow,
	Wall, 
	Sym
};

// Supersonic Inflow
class t_BCDataInflow :public t_BCDataFace {
public:

	static const std::string bc_kind;

	t_BCDataInflow() = delete;
	t_BCDataInflow(const std::string& sect) :t_BCDataFace(sect) { default_settings(); }
	const std::string& getBCKindName() const { return bc_kind; }
	// implement TPugin
	std::string get_name() const { return bc_kind + "/" + nameOfFldSection; }
	std::string get_description() const { return std::string("bc inflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// this is supersonic outflow, i.e. extrapolation
class t_BCDataOutFlow :public t_BCDataFace {

public:

	static const std::string bc_kind;

	t_BCDataOutFlow() = delete;
	t_BCDataOutFlow(const std::string& sect) :t_BCDataFace(sect) { default_settings(); }
	const std::string& getBCKindName() const { return bc_kind; }
	// implement TPugin
	std::string get_name() const { return bc_kind + "/" + nameOfFldSection; }
	std::string get_description() const { return std::string("bc outflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// Euler wall
class t_BCDataEulerWall :public t_BCDataFace {

	double TWallDim;

public:

	static const std::string bc_kind;

	t_BCDataEulerWall() = delete;
	t_BCDataEulerWall(const std::string& sect) :t_BCDataFace(sect) { default_settings(); }
	const std::string& getBCKindName() const { return bc_kind; }
	// implement TPugin
	std::string get_name() const { return bc_kind + "/" + nameOfFldSection; }
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
	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// Symmetry bc
class t_BCDataSym :public t_BCDataFace {
public:

	static const std::string bc_kind;

	t_BCDataSym() = delete;
	t_BCDataSym(const std::string& sect) :t_BCDataFace(sect) { default_settings(); }
	const std::string& getBCKindName() const { return bc_kind; }
	// implement TPugin
	std::string get_name() const { return bc_kind + "/" + nameOfFldSection; }
	std::string get_description() const { return std::string("bc sym"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

/**
 * List of bc sets for euler eqs
 // read bc_list section from config dynamically
 // so use iniFile capabilities directly instead of plugin params
 */

class t_BCListEuler : public t_BCList {

public:
	// implement TPugin
	t_BCListEuler() { default_settings(); }
	std::string get_name() const { return "bc_list"; };
	std::string get_description() const { return std::string("bc list"); };
	void default_settings() {	};
	// implement BCList
	std::string getSupportedBCsStr();
	t_FaceBCID getBCID(std::string sectionname);
	bool getBCKindBySectName(const std::string& name, t_BCKindEuler& kind);
	void addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data);
	bool has(std::string sectionname);
	void init(std::string& ini_data, const std::string& spec) {

		t_BCList::init(ini_data, spec);

	}

};

extern t_BCListEuler G_BCListEuler;

