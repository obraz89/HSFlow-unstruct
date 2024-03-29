#include "settings.h"

#include "logging.h"

#include <cassert>

#if defined(_WINDOWS)
#include <direct.h>  // mkdir
#define stat  _stat64
#define getcwd  _getcwd
#else
#include <unistd.h>
#endif

const char* g_CASE_RESULTS_DIR = "solution/";   // ending '/' is required!
const char* g_CASE_SETTINGS_DIR = "settings/";

static const char MAIN_INI[] = "main.ini";

t_Settings g_Settings;

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

bool hs_file_exists(const std::string& fn)
{
	struct stat buf;
	if (stat(fn.c_str(), &buf) == 0)
		return bool(buf.st_mode & S_IFREG);

	return false;
}

bool hs_dir_exists(const std::string& dn)
{
	struct stat buf;
	if (stat(dn.c_str(), &buf) == 0)
		return bool(buf.st_mode & S_IFDIR);

	return false;
}

bool hs_dir_create(const std::string& dn)
{
	if (hs_dir_exists(dn))
		return true;

	int rc = 0;

#if defined(_WINDOWS)
	rc = _mkdir(dn.c_str());
#else
	mode_t nMode = 0777; // UNIX style permissions, will be reduced by umask
	rc = mkdir(dn.c_str(), nMode);
#endif

	return rc == 0;
}

bool load_settings() {

	const std::string sIni = std::string() + g_CASE_SETTINGS_DIR + MAIN_INI;

	short ok = 1;

	// Create config & solution directories
	if (true)
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
			std::string& fn = g_Settings.strGridFnUstr;
			fn = iniAD.read_string("GridFileUnstruct", "");
			ok &= CheckFieldFile(fn);  // don't return here, allow to init the whole config
		}

		{
			std::string& fn = g_Settings.strFieldFnUstr;
			fn = iniAD.read_string("FieldFileUnstruct", "");
			ok &= CheckFieldFile(fn);  // don't return here, allow to init the whole config
		}

		{
			std::string& fn = g_Settings.strGridFnStrc;
			fn = iniAD.read_string("GridFileStruct", "");
			ok &= CheckFieldFile(fn);  // don't return here, allow to init the whole config
		}

	}

	ini_data = iniAD.get_ini().encode();

	iniAD.get_ini().decode(ini_data);

	iniAD.get_ini().save(sIni);

	return ok;

}