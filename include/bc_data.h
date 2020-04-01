#pragma once

#include <map>

#include "PluginBase.h"

#include "common_data.h"

#define BC_INFLOW_STR "bc_inflow"
#define BC_OUTFLOW_STR "bc_outflow"
#define BC_EULER_WALL_STR "bc_euler_wall"
#define BC_SYM_STR "bc_sym"

// base class for information of face-type boundary conditions
class t_BCDataFace : public TPlugin{
protected:
	std::string name;
public:
	t_BCDataFace(const std::string& n) :name(n) {}
	std::string get_name() const { return name; };
	virtual void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) = 0;

};

// Inflow, VelInf should be non-dim unit vector
// zero angle of attack case: (1,0,0)
class t_BCDataInflow :public t_BCDataFace {

	t_Vec3 VelInf;

public:
	t_BCDataInflow() :t_BCDataFace(BC_INFLOW_STR) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc inflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// this is supersonic outflow, i.e. extrapolation
class t_BCDataOutFlow :public t_BCDataFace {

public:
	t_BCDataOutFlow() :t_BCDataFace(BC_OUTFLOW_STR) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc outflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// Euler wall
class t_BCDataEulerWall :public t_BCDataFace {

	double TWallDim;

public:
	t_BCDataEulerWall() :t_BCDataFace(BC_EULER_WALL_STR) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc euler wall"); };
	void default_settings() {
		TPluginParamsGroup g("", "gas-dynamic functions values on the wall");
		g.add("Tw_K", 300.0, "temperature, dimensional in K");
		_mapParamsGrps.emplace(g.get_name(), g);
	};
	virtual void init(std::string& ini_data, const std::string& spec) {

		TPlugin::init(ini_data, spec);

		const TPluginParamsGroup& g = get_settings_grp("");

		TWallDim = g.get_real_param("Tw_K");

	};

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// Symmetry bc
class t_BCDataSym :public t_BCDataFace {
public:
	t_BCDataSym() :t_BCDataFace(BC_SYM_STR) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("bc sym"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// This is identical to BC set in  HSFlow
// TODO: if necessary, can be a dynamic list
struct t_BCList {
	t_BCDataInflow bc_inflow;
	t_BCDataOutFlow bc_outflow;
	t_BCDataEulerWall bc_euler_wall;
	t_BCDataSym bc_sym;
	// to iterate
	std::vector<t_BCDataFace*> pBCs;
	t_BCList() :bc_inflow(), bc_outflow(), bc_euler_wall(), bc_sym(), pBCs() {
		pBCs.push_back(&bc_inflow);
		pBCs.push_back(&bc_outflow);
		pBCs.push_back(&bc_euler_wall);
		pBCs.push_back(&bc_sym);
	}
	const t_BCDataFace* getBCByKind(t_FaceBCKind kind) {

		switch (kind) {
		case t_FaceBCKind::Inflow:
			return &bc_inflow;
			break;
		case t_FaceBCKind::Outflow:
			return &bc_outflow;
			break;
		case t_FaceBCKind::Wall:
			return &bc_euler_wall;
			break;
		case t_FaceBCKind::Sym:
			return &bc_sym;
			break;
		default:
			hsLogMessage("Error:t_BCList: bc handler was not loaded for this bc");
			break;
		}

	}
	void init(std::string& ini_data, const std::string& spec) {

		for (int i = 0; i < pBCs.size(); i++)
			pBCs[i]->init(ini_data, spec);

	}
};

extern t_BCList G_BCList;
