#pragma once

#include "IniFile.hpp"

#include "logging.h"

#include <cassert>

//
// Constants
//
extern const char* g_CASE_SETTINGS_DIR;
extern const char* g_CASE_RESULTS_DIR;

// enum int with str code
// 0 <=> Str1
// 1 <=> Str2
// etc
struct t_EnumStr {

	int Val;

	std::vector<std::string> ValsStr;

	virtual void initValsStr() = 0;

	void set(std::string str) {
		for (int i = 0; i < ValsStr.size(); i++) {
			if (ValsStr[i].compare(str) == 0) {
				Val = i;
				return;
			}
		}
		hsLogError("EnumStr: unknown option %s", str.c_str());
	}

	std::string getOptionsStr() {
		std::string str;
		for (auto elem : ValsStr) str += elem + ",";
		return str;
	}

	std::string toStr() {
		if (Val < ValsStr.size()) return ValsStr[Val];
		hsLogError("t_EnumStr:wrong val");
		return "";
	}

	bool operator==(int v) { return Val == v; }

	std::string defaultValStr() const { 
		return ValsStr.size()>0 ? ValsStr[0] : ""; }

};

struct t_NonDimType : public t_EnumStr {
	static const int FreeStreamVelo = 0;
	static const int FreeStreamSoundSpeed = 1;
	void initValsStr() {
		ValsStr.push_back("FreeStreamVelo");
		ValsStr.push_back("FreeStreamSoundSpeed");
	}

	t_NonDimType() { initValsStr(); }
};

struct t_EnumInitFieldCustom : public t_EnumStr {
	static const int No = 0;
	static const int Eu1d = 1;
	void initValsStr() {
		ValsStr.push_back("No");
		ValsStr.push_back("Eu1d");
	}

	t_EnumInitFieldCustom() { initValsStr(); }
};

struct TgenericSettings
{
	//--- Init ---
	std::string strInitFieldFN; // init field filename
//	std::string strPrevFieldFN; // field at previous time step, approximating if empty
	std::string strGridFN;      // grid file name
//	bool allowFrozenZones;      // enable frozen zones if specified in the grid
//	float init_ltzero_filter; // filter out non-physical negative values replacing by `init_ltzero_filter * U_inf`, <=0 - if not needed

	// if field to be initialized by custom procedure (e.g. euler 1d)
	t_EnumInitFieldCustom initFieldCustom;
	//--- Time ---
	double timeStart;  // initial time, if <0 - take from init-file
	int numTimeSteps; // number of time steps to do
	double convergenceResidual; // stop solver at this residual

	// euler, ns, rans etc
	std::string strCase;

	// non-dim type
	t_NonDimType nonDimType;

	// TODO: scheme (domain) options
	double CFL;	// Courant number
	std::string strScheme;
	std::string strRiemannSolver;

	// Automatic increasing of time step size:
//	struct {
//		double factor;       // dt *= factor
//		int nSteps_interval; // check after this number of steps
//		double limit;        // increase timestep up to this value
//		double minResidual;  // increase timestep only if 0-th NES residual is lower
//	} dt_inc;

	//--- Iterations ---
//	double epsOut;   //Newton iteration residual
//	int iterOut;     //max num of Newton iterations
//	double epsIn;    //internal iterations residual (GMRES)
//	int iterIn;      //max num of internal iterations (GMRES)
//	int itsFalseJac; // max num of iterations with False Jacoby matrix
//	double startNwtM;   //init value for tau param in modified Newton method
//	double qmax;		//relative decrement of residual at Newton step -- indication to recal Jac matrix

	//--- Saving ---
	std::string strSaveFieldFNTemplate;  // printf-style template for generating file name using time
	std::string strSaveGridFN;           // name of separate grid file, if empty grid is embedded

	int timeSteps2Write;            // time steps count to save full field
	//bool isDoublePrecision2Write;   // precision to write full field

//	int slices_timeSteps2Write;    // time steps count to write field slices (<=0 if not saving)
//	std::vector<TSliceZneDef> lstSlices;
	//bool slices_isDoublePrecision2Write;  // precision to write slices

//	int wall_timeSteps2Write;          // time steps count to write wall data (<=0 if not saving)
//	std::set<std::string> setWallFuncs;

//  bool wall_isDoublePrecision2Write;   // precision to write wall
};

extern TgenericSettings g_genOpts;


//-----------------------------------------------------------------------------

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
		catch (std::domain_error& e) {
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
		catch (std::domain_error& e) {
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
		catch (std::domain_error& e) {
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
		catch (std::domain_error& e) {
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
		catch (std::domain_error& e) {
			hsLogWarning("%s", e.what());
		}

		sect[key] = b0;   _updated = true;
		return b0;
	}
};

bool load_settings();

void load_case(std::string& ini_data);

void load_case_euler(std::string& ini_data);

void load_case_ns(std::string& ini_data);