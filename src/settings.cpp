#include "settings.h"

#include "common_data.h"
#include "common_procs.h"

#include "flow_params.h"

#include "bc_common.h"

#include "mesh.h"

#include <stdlib.h> // free()
#include <string.h> // strtok(), strdup()
#include <cassert>

const char* g_CASE_RESULTS_DIR = "solution/";   // ending '/' is required!
const char* g_CASE_SETTINGS_DIR = "settings/";
//-----------------------------------------------------------------------------

static const char MAIN_INI[] = "main.ini";

TgenericSettings g_genOpts;

/**
 * Wrapper class for IniFile
 * Writes default value if key doesn't exist on reading
 */
class TIniAutoDefaults
{
	ini::IniFile _ini;
	std::string _file_name;

	ini::IniSection* _section;  // active section
	std::string _section_name;

	bool _updated;
	TIniAutoDefaults();

public:
	TIniAutoDefaults(const std::string& fn)
		: _ini(fn), _file_name(fn), _section(nullptr), _section_name(), _updated(false) {
		;
	}

	ini::IniFile& get_ini() {
		return _ini;
	}

	void save_if_updated() {
		if (_updated)
			_ini.save(_file_name);
	}

	void set_section(const char* s) {
		_section_name = s;
		_section = &(_ini[s]);
	}

	/**
	 * Read the key or create a new one with defaut value
	 *
	 * @param[in] key - key name in the active section
	 * @param[in] s0 - default value in case the key doesn't not exist
	 * @return    key value
	 */
	std::string read_string(const std::string& key, const std::string& s0) {
		assert(_section);
		ini::IniSection& sect = *_section;
		if (sect.has(key))
			return sect[key].asString();

		sect[key] = s0;   _updated = true;
		return s0;
	}

	int read_int(const std::string& key, int i0) {
		assert(_section);
		ini::IniSection& sect = *_section;
		try {
			if (sect.has(key))
				return sect[key].asInt();
		}
		catch (std::domain_error & e) {
			hsLogWarning("%s", e.what());
		}

		sect[key] = i0;   _updated = true;
		return i0;
	}

	/**
	 * Reread the key, if updated assign the provided variable
	 *
	 * @param[in]     key  - key name in the active section
	 * @param[in/out] prevValue - saved previous value of the key
	 * @return        true if value was updated and false otherwise
	 */
	bool reread_int(const std::string& key, int& prevValue) {
		int val;
		try {
			if (!_section->has(key))
				return false;
			val = (*_section)[key].asInt();
		}
		catch (std::domain_error & e) {
			return false;
		}

		if (val == prevValue)
			return false;

		hsLogMessage("  %s/%s: %d -> %d", _section_name.c_str(), key.c_str(), prevValue, val);
		prevValue = val;

		return true;
	}

	double read_float(const std::string& key, double f0) {
		assert(_section);
		ini::IniSection& sect = *_section;
		try {
			if (sect.has(key))
				return sect[key].asDouble();
		}
		catch (std::domain_error & e) {
			hsLogWarning("%s", e.what());
		}

		sect[key] = f0;   _updated = true;
		return f0;
	}

	/**
	 * Reread the key, if updated assign the provided variable
	 *
	 * @param[in]     key  - key name in the active section
	 * @param[in/out] prevValue - saved previous value of the key
	 * @return        true if value was updated and false otherwise
	 */
	bool reread_float(const std::string& key, double& prevValue) {
		double val;
		try {
			if (!_section->has(key))
				return false;
			val = (*_section)[key].asDouble();
		}
		catch (std::domain_error & e) {
			return false;
		}

		if (fabs(val - prevValue) < 1e-15)
			return false;

		hsLogMessage("  %s/%s: %g -> %g", _section_name.c_str(), key.c_str(), prevValue, val);
		prevValue = val;

		return true;
	}

	bool read_bool(const std::string& key, bool b0) {
		assert(_section);
		ini::IniSection& sect = *_section;
		try {
			if (sect.has(key))
				return sect[key].asBool();
		}
		catch (std::domain_error & e) {
			hsLogWarning("%s", e.what());
		}

		sect[key] = b0;   _updated = true;
		return b0;
	}
};


bool load_settings() {

	const std::string sIni = std::string() + g_CASE_SETTINGS_DIR + MAIN_INI;

	G_State.mpiNProcs = 1;

	G_State.mpiRank = 0;

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

		g_genOpts.strGridFN = iniAD.read_string("gridFN", "grid.cgns");
		g_genOpts.timeStart = iniAD.read_float("timeStart", -1);

		// time steps & flow saving
		g_genOpts.numTimeSteps = iniAD.read_int("NumTimeSteps", 100000);
		g_genOpts.timeSteps2Write = iniAD.read_int("TimeSteps2Write", 500);

	}

	iniAD.set_section("scheme");
	{
		g_genOpts.CFL = iniAD.read_float("CFL", 0.1);
		g_genOpts.strScheme = iniAD.read_string("domain", "euler_1st");
	}



	// TODO: just do this in MPI case
	//TPlugin::load_settings(fn, ini_data);

	ini_data = iniAD.get_ini().encode();

	G_FreeStreamParams.init(ini_data, "");

	G_pBCList->init(ini_data, "");

	// for now rewriting ini every time (to get defaults for the first time) 
	// TODO: replace by TPlugin::save_settings(fn, ini_data); 
	// when MPI is ok
	iniAD.get_ini().decode(ini_data);

	iniAD.get_ini().save(sIni);

	return true;

};