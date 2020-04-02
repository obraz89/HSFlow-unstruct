#pragma once

#include "PluginBase.h"

#include "common_data.h"

#include <vector>

// base class for information of face-type boundary conditions
class t_BCDataFace : public TPlugin{
protected:
	std::string nameOfFldSection;
public:
	t_BCDataFace() = delete;
	t_BCDataFace(const std::string& sect):nameOfFldSection(sect) {}
	virtual void yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) = 0;
	virtual ~t_BCDataFace() {}

};

// Supersonic Inflow
class t_BCDataInflow :public t_BCDataFace {
public:

	static const std::string bc_kind;

	t_BCDataInflow() = delete;
	t_BCDataInflow(const std::string& sect):t_BCDataFace(sect){ default_settings(); }

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
	t_BCDataOutFlow(const std::string& sect):t_BCDataFace(sect){ default_settings(); }

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
	t_BCDataSym(const std::string& sect):t_BCDataFace(sect){ default_settings(); }
	// implement TPugin
	std::string get_name() const { return bc_kind + "/" + nameOfFldSection; }
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
			g.add(&(bc_sets[i][0]), t_BCDataInflow::bc_kind, "bc set");
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

				t_BCDataFace* pBC = nullptr;

				if (bc_kind_name.compare(t_BCDataInflow::bc_kind) == 0)
					pBC = new t_BCDataInflow(bc_set_name);

				if (bc_kind_name.compare(t_BCDataOutFlow::bc_kind) == 0)
					pBC = new t_BCDataOutFlow(bc_set_name);

				if (bc_kind_name.compare(t_BCDataEulerWall::bc_kind) == 0)
					pBC = new t_BCDataEulerWall(bc_set_name);

				if (bc_kind_name.compare(t_BCDataSym::bc_kind) == 0)
					pBC = new t_BCDataSym(bc_set_name);

				if (pBC==nullptr)
					hsLogMessage("Error:t_BCList:Unknown bc type!");
				else {
					hsLogMessage("Initializing bc set: %s", &(pBC->get_name()[0]));
					pBC->init(ini_data, "");
					_pBCs.emplace(std::make_pair(bc_set_name, pBC));
				}


			}

	}
	std::string getSupportedBCsStr();

	~t_BCList() {
		std::map<std::string, t_BCDataFace*>::iterator it;

		for (auto elem : _pBCs) delete elem.second;

	}

};



	//void init_boco(const std::string& patchFamily);

extern t_BCList G_BCList;