#pragma once

#include "PluginBase.h"

#include "common_data.h"

#include <vector>

#define BC_INFLOW_STR "bc_inflow"
#define BC_OUTFLOW_STR "bc_outflow"
#define BC_EULER_WALL_STR "bc_euler_wall"
#define BC_SYM_STR "bc_sym"

// base class for information of face-type boundary conditions
class t_BCDataFace : public TPlugin{
protected:
	std::string nameOfFldSection;
public:
	void setFldSectionName(const std::string& section_name) { nameOfFldSection = section_name; };
	std::string getFldSectionName() { return nameOfFldSection; }
	virtual void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) = 0;

};

// Supersonic Inflow
class t_BCDataInflow :public t_BCDataFace {

public:
	t_BCDataInflow(){ default_settings(); }
	// implement TPugin
	std::string get_name() const { return BC_INFLOW_STR; };
	std::string get_description() const { return std::string("bc inflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// this is supersonic outflow, i.e. extrapolation
class t_BCDataOutFlow :public t_BCDataFace {

public:
	t_BCDataOutFlow(){ default_settings(); }
	// implement TPugin
	std::string get_name() const { return BC_OUTFLOW_STR; };
	std::string get_description() const { return std::string("bc outflow"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

// Euler wall
class t_BCDataEulerWall :public t_BCDataFace {

	double TWallDim;

public:
	// implement TPugin
	std::string get_name() const { return BC_EULER_WALL_STR; };
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
	t_BCDataSym(){ default_settings(); }
	// implement TPugin
	std::string get_name() const { return BC_SYM_STR; };
	std::string get_description() const { return std::string("bc sym"); };
	void default_settings() {};
	void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };

	void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs);
};

std::vector<std::string> tokenize_str(std::string str);

/**
 * Boundary condition sub-plugins manager,
 * handles plugins references to each-other
 */
class t_BCList : public TPlugin{

	std::map<std::string, t_BCDataFace*>  _pBCs;

public:
	// implement TPugin
	t_BCList() { default_settings(); }
	std::string get_name() const { return "bc_list"; };
	std::string get_description() const { return std::string("bc list"); };
	void default_settings() {

		TPluginParamsGroup g("", "list of bcs");

		g.add("Supported_bc_names", getSupportedBCsStr(), "supported kinds of bcs");
		g.add("bc_list", "inflow, outflow, wall, sym", "list of bcs");

		std::vector<std::string> bc_sets = tokenize_str(g.get_string_param("bc_list"));

		for (int i = 0; i < bc_sets.size(); i++) {
			g.add(&(bc_sets[i][0]), BC_INFLOW_STR, "bc set");
		}

		_mapParamsGrps.emplace(g.get_name(), g);


	};
	void init(std::string& ini_data, const std::string& spec) { 

		TPlugin::init(ini_data, spec); 

		const TPluginParamsGroup& g = get_settings_grp("");

		std::vector<std::string> bc_sets = tokenize_str(g.get_string_param("bc_list"));

			for (int i = 0; i < bc_sets.size(); i++) {
				std::string bc_set_name = bc_sets[i];
				std::string bc_kind_name = g.get_string_param(&(bc_sets[i][0]));
				hsLogMessage("Addind BC with sec_name %s : type is %s", 
					&bc_set_name[0], &bc_kind_name[0]);
			}

	};

	std::string getSupportedBCsStr();

	//void init_boco(const std::string& patchFamily);
};

extern t_BCList G_BCList;