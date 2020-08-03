#pragma once

#include "IniFile.hpp"

//
// Constants
//
extern const char* g_CASE_SETTINGS_DIR;
extern const char* g_CASE_RESULTS_DIR;

struct t_Settings {

	std::string strGridFnUstr;
	std::string strGridFnStrc;

	std::string strFieldFnUstr;
};

extern t_Settings g_Settings;

bool load_settings();
