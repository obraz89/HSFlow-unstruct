#include "settings.h"

#include "common_data.h"

#include "flow_params.h"

#include "bc_data.h"


void load_settings(const std::string& fn) {

	G_Domain.nDim = 3;

	G_State.mpiNProcs = 1;

	G_State.mpiRank = 0;

	std::string ini_data;

	// TODO: just do this in MPI case
	//TPlugin::load_settings(fn, ini_data);

	ini::IniFile ini;
	try {
		ini.load(fn);
	}
	catch (std::logic_error & e) {
		hsLogError("Can't load settings from '%s' (%s)", fn.c_str(), e.what());
	}

	ini_data = ini.encode();

	G_FreeStreamParams.init(ini_data, "");

	G_BCList.init(ini_data, "");

	// for now rewriting ini every time (to get defaults for the first time) 
	// TODO: replace by TPlugin::save_settings(fn, ini_data); 
	// when MPI is ok
	ini.decode(ini_data);

	ini.save(fn);

};