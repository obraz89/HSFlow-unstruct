#include "settings.h"

#include "common_data.h"
#include "common_procs.h"

#include "flow_params.h"

#include "bc_common.h"

#include "dom_base.h"

#include "flow_model_perfect_gas.h"

#include <stdlib.h> // free()
#include <string.h> // strtok(), strdup()

const char* g_CASE_RESULTS_DIR = "solution/";   // ending '/' is required!
const char* g_CASE_SETTINGS_DIR = "settings/";
//-----------------------------------------------------------------------------

static const char MAIN_INI[] = "main.ini";

TgenericSettings g_genOpts;


bool load_settings() {

	const std::string sIni = std::string() + g_CASE_SETTINGS_DIR + MAIN_INI;

	short ok = 1;

	// Create config & solution directories
	if (G_State.mpiRank == 0)
	{
		ok &= hs_dir_create(g_CASE_SETTINGS_DIR) ? 1 : 0;
		ok &= hs_dir_create(g_CASE_RESULTS_DIR) ? 1 : 0;

		if (!ok)
			hsLogError("Can't create '%s' and/or '%s' case dirs",
				g_CASE_SETTINGS_DIR, g_CASE_RESULTS_DIR);
	}

	TIniAutoDefaults iniAD(sIni);

	std::string ini_data;

	iniAD.set_section("init");
	{
		struct {
			short operator()(std::string& fn) {
				if (fn.empty())  return 1;
				if (hs_file_exists(fn))  return 1;

				const std::string fn2 = g_CASE_RESULTS_DIR + fn;
				if (hs_file_exists(fn2)) { fn = fn2;  return 1; }

				hsLogError("None of '%s', '%s' field files were found", fn.c_str(), fn2.c_str());
				return 0;
			}
		}  CheckFieldFile;

		{
			std::string& fn = g_genOpts.strInitFieldFN;
			fn = iniAD.read_string("initFieldFN", "");
			ok &= CheckFieldFile(fn);  // don't return here, allow to init the whole config
		}

		// initialize field with custom procedure
		{
			iniAD.read_string("initFieldCustom_Options",
				g_genOpts.initFieldCustom.getOptionsStr());

			g_genOpts.initFieldCustom.set(iniAD.read_string("InitFieldCustom",
				g_genOpts.initFieldCustom.defaultValStr()));
		}


		g_genOpts.strGridFN = iniAD.read_string("gridFN", "grid.cgns");
		g_genOpts.timeStart = iniAD.read_float("timeStart", -1);

		// time steps & flow saving
		g_genOpts.numTimeSteps = iniAD.read_int("NumTimeSteps", 100000);
		g_genOpts.timeSteps2Write = iniAD.read_int("TimeSteps2Write", 500);

		g_genOpts.strCase = iniAD.read_string("case", "euler");

		iniAD.read_string("nonDimOptions", g_genOpts.nonDimType.getOptionsStr());

		g_genOpts.nonDimType.set(iniAD.read_string("nonDim", 
			g_genOpts.nonDimType.defaultValStr()));

	}

	iniAD.set_section("scheme");
	{
		g_genOpts.CFL = iniAD.read_float("CFL", 0.1);
		g_genOpts.strScheme = iniAD.read_string("domain", "euler_1st");
		g_genOpts.strRiemannSolver = iniAD.read_string("RiemannSolver", "Roe");
	}

	// TODO: just do this in MPI case
	//TPlugin::load_settings(fn, ini_data);
	ini_data = iniAD.get_ini().encode();

	load_case(ini_data);

	iniAD.set_section("flow_model");

	{
		G_FlowModelParams.Gamma = 1.4;

		G_FlowModelParams.Pr = 0.72;

		iniAD.read_string("ViscTypeOptions", G_FlowModelParams.ViscType.getOptionsStr());

		G_FlowModelParams.ViscType.set(iniAD.read_string("ViscType", 
			G_FlowModelParams.ViscType.defaultValStr()));

	}

	// for now rewriting ini every time (to get defaults for the first time) 
	// TODO: replace by TPlugin::save_settings(fn, ini_data); 
	// when MPI is ok
	ini_data = iniAD.get_ini().encode();
	iniAD.get_ini().decode(ini_data);
	iniAD.get_ini().save(sIni);

	return true;

};

void load_case(std::string& ini_data) {

	if (g_genOpts.strCase.compare("euler") == 0) {

		load_case_euler(ini_data);
		return;

	}

	if (g_genOpts.strCase.compare("ns") == 0) {

		load_case_ns(ini_data);
		return;

	}

	hsLogMessage("Error: unsupported case option, supported options are: euler, ns");

}