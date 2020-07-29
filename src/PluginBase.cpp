///////////////////////////////////////////////////////////////////////////////
// Name:        PluginBase.cpp
// Purpose:     Default behaviour of base classes for plugins
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include <mpi.h>

// Dynamic library support
#if defined(_WINDOWS)
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <stdlib.h>  // getenv()

#include "IniFile.hpp"

#include "PluginBase.h"

#include "common_data.h"

static const char HSFLOW_ENV_VAR[] = "HSFLOW";
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class TPlugin
// NB: Keep definitions in *.cpp unit to prevent leaking into external plugins
///////////////////////////////////////////////////////////////////////////////

TPlugin::TPlugin()
{
	default_settings();
}

TPlugin::~TPlugin() = default;
//-----------------------------------------------------------------------------

TPluginParamsGroup& TPlugin::get_settings_grp(const char* szGrpName)
{
	try {
		return _mapParamsGrps.at(szGrpName);
	}
	catch (const std::out_of_range&) {
		hsTHROW("%s: Settings group '%s' doesn't exist", get_name().c_str(), szGrpName);
	}
}

const TPluginParamsGroup& TPlugin::get_settings_grp(const char* szGrpName) const
{
	try {
		return _mapParamsGrps.at(szGrpName);
	}
	catch (const std::out_of_range&) {
		hsTHROW("%s: Settings group '%s' doesn't exist", get_name().c_str(), szGrpName);
	}
}
//-----------------------------------------------------------------------------


/**
 * Default initializing of a plugin
 * WARNING: NOT COLLECTIVE function
 *
 * @param[in] ini_data - config inifile data
 * @param[in] spec - specialization of the plugin, e.g. face name it belongs to
 */
void TPlugin::init(std::string & ini_data, const std::string & spec)
{
	_spec = spec;
	parse_settings(ini_data);
}
//-----------------------------------------------------------------------------


/**
 * Loads settings from inifile into a string & broadcasts it to every MPI rank
 * COLLECTIVE operation (should be called by every MPI rank)
 *
 * @param[in]  fn - INI-file name to load settings from
 * @param[out] ini_data - memory buffer to store setting data
 *
 * @return true if succeded, false otherwise
 */
bool TPlugin::load_settings(const std::string & fn, std::string & ini_data)
{
	if (G_State.mpiRank == 0)
	{
		ini::IniFile ini;
		try {
			ini.load(fn);
		}
		catch (std::logic_error & e) {
			hsLogError("Can't load settings from '%s' (%s)", fn.c_str(), e.what());
		}

		ini_data = ini.encode();
	}

	int data_len = ini_data.length();

	MPI_Bcast(&data_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (data_len <= 0)
		return false;

	ini_data.resize(data_len);
	MPI_Bcast(&ini_data[0], data_len, MPI_CHAR, 0, MPI_COMM_WORLD);
	return true;
}
//-----------------------------------------------------------------------------


/**
 * Saves settings gathered from every MPI rank into inifile @ root
 * COLLECTIVE operation (should be called by every MPI rank)
 *
 * @param[in]  fn - INI-file name to save settings into
 * @param[out] ini_data - memory buffer havind data to save
 */
void TPlugin::save_settings(const std::string & fn, const std::string & ini_data)
{
	const int tag_len = 'i' + 'n' + 'i' + 'L';
	const int tag_data = 'i' + 'n' + 'i' + 'D';

	if (G_State.mpiRank == 0)
	{
		// Original inifile at root rank to merge new data into
		ini::IniFile ini0(fn);

		// Merge in new sections and/or fields
		auto merge = [&ini0](const std::string& data)
		{
			bool is_updated = false;

			ini::IniFile ini;   ini.decode(data);
			for (auto& sect : ini) {
				if (!ini0.has(sect.first)) {
					// Copy new section to root ini
					ini0[sect.first] = sect.second;
					is_updated = true;
				}
				else {
					ini::IniSection& sect0 = ini0[sect.first];
					for (const auto& field : sect.second) {
						// Copy new fields to existing section
						if (!sect0.has(field.first)) {
							sect0[field.first] = field.second;
							is_updated = true;
						}
					}
				}
			}
			return is_updated;
		};

		// Merge in updated data at root rank
		bool is_updated = merge(ini_data);

		// Recieve data from other ranks & merge in
		for (int p = 1; p < G_State.mpiNProcs; ++p)
		{
			int data_len = 0;
			MPI_Recv(&data_len, 1, MPI_INT, p, tag_len, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			std::string data(data_len, '\0');
			MPI_Recv(&data[0], data_len, MPI_CHAR, p, tag_data, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			is_updated |= merge(data);
		}

		if (is_updated)
			ini0.save(fn);
	}
	else
	{
		// Send data to root MPI rank
		const int data_len = ini_data.length();
		MPI_Ssend(const_cast<int*>(&data_len), 1, MPI_INT, 0, tag_len, MPI_COMM_WORLD);

		MPI_Ssend(const_cast<char*>(&ini_data[0]), data_len, MPI_CHAR, 0, tag_data, MPI_COMM_WORLD);
	}
}
//-----------------------------------------------------------------------------


/**
 * Parse plugin settings from ini-file data
 * WARNING: Function is NOT collective (plugins may call it different number of times)
 *
 * @param[in/out] ini_data - encoded ini-file data, may be empty to be filled
 */
void TPlugin::parse_settings(std::string & ini_data)
{
	if (_mapParamsGrps.empty())  // plugin has no parameters
		return;

	ini::IniFile ini;   ini.decode(ini_data);
	bool isIniUpdated = false;

	std::string section = get_name();

	// Add ini-section for the plugin if not exist
	if (!ini.has(section)) {
		ini[section].comment = get_description();
		isIniUpdated = true;
	}

	if (!_spec.empty())  section += "/" + _spec;
	for (auto& grp : _mapParamsGrps)
	{
		const std::string& group = grp.first;
		const std::string subsection = section + (group.empty() ? "" : "/" + group);

		for (auto& prm : grp.second.mapParams)
		{
			const std::string& key = prm.first;
			TPluginParam& pp = prm.second;
			const std::string& val0 = pp.get_raw_value();

			ini::IniSection& iniSect = ini[subsection];
			if (!iniSect.has(key)) {
				// Key not found -> add it to ini-file with default value
				hsLogWarning("Settings file has no [%s]/%s. Assigning '%s'",
					subsection.c_str(), key.c_str(), val0.c_str());

				//const std::string& comment = pp.description;
				iniSect[key] = val0;
				isIniUpdated = true;
			}
			else {
				// Key found -> store its value into TPluginParam
				if (!pp.set_raw_value(iniSect[key].asString())) {
					hsLogWarning("Wrong value for [%s]/%s. Using '%s'", subsection.c_str(), key.c_str(), val0.c_str());
				}
			}
		} // loop params
	} // loop groups

	if (isIniUpdated)
		ini_data = ini.encode();
}
//-----------------------------------------------------------------------------

// Explicitly instantiate the class variants to be used
template class OrderedMap<std::string, TPluginParam>;
template class OrderedMap<std::string, TPluginParamsGroup>;
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class TPluginParamsGroup
// NB: Keep definitions in *.cpp unit to prevent leaking into external plugins
///////////////////////////////////////////////////////////////////////////////

TPluginParamsGroup::TPluginParamsGroup() = default;

TPluginParamsGroup::TPluginParamsGroup(const char* aName, const char aDescr[])
	: m_name(aName), m_description(aDescr)
{
	;
}

TPluginParamsGroup::~TPluginParamsGroup() = default;
//-----------------------------------------------------------------------------

const TPluginParam&
TPluginParamsGroup::get_raw_param(const std::string & parName) const
{
	try {
		return mapParams.at(parName);
	}
	catch (const std::out_of_range&) {
		hsTHROW("Unable to find parameter '%s' in the group '%s'",
			parName.c_str(), m_name.c_str());
	}
}
//-----------------------------------------------------------------------------

double TPluginParamsGroup::get_real_param(const char* szName) const
{
	const TPluginParam& param = get_raw_param(szName);

	double val;
	if (!param.get_value(val))
		hsTHROW(
			"Parameter '%s' in the group '%s' has value '%s', which is not a real number.",
			szName, m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

int TPluginParamsGroup::get_int_param(const char* szName) const
{
	const TPluginParam& param = get_raw_param(szName);

	int val;
	if (!param.get_value(val))
		hsTHROW(
			"Parameter '%s' in the group '%s' has value '%s', which is not an integer.",
			szName, m_name.c_str(), param.get_raw_value().c_str()
		);

	return val;
}
//-----------------------------------------------------------------------------

const std::string& TPluginParamsGroup::get_string_param(const char* szName) const
{
	const TPluginParam& par = get_raw_param(szName);

	const std::string* val = nullptr;
	if (!par.get_value(val))
		hsTHROW(
			"Parameter '%s' in the group '%s' has value '%s', which is not of string type.",
			szName, m_name.c_str(), par.get_raw_value().c_str()
		);

	return *val;
}

std::vector<std::string> TPluginParamsGroup::get_list_of_names() const {

	return mapParams.get_list_of_keys();

};
//-----------------------------------------------------------------------------


///////////////////////////////////////////////////////////////////////////////
// class TPluginParam
// NB: Keep definitions in *.cpp unit to prevent leaking into external plugins
///////////////////////////////////////////////////////////////////////////////

bool TPluginParam::set_raw_value(const std::string & val)
{
	// Check the value
	switch (type)
	{
	case ptInt:
		int n;   if (!(std::istringstream(val) >> n)) {
			hsLogWarning("Provided value '%s' is not of required type 'int'", val.c_str());
			return false;
		}
		break;

	case ptDouble:
		double x;   if (!(std::istringstream(val) >> x)) {
			hsLogWarning("Provided value '%s' is not of required type 'double'", val.c_str());
			return false;
		}
		break;

	case ptString:
		break; // fallthrough
	}

	_value = val;  //TODO: _value.Trim();
	return true;
}
//-----------------------------------------------------------------------------


bool TPluginParam::set_value(int val)
{
	if (type != ptInt)  return false;

	std::ostringstream ss;  ss << val;
	if (ss.fail()) return false;

	_value = ss.str();
	return true;
}

bool TPluginParam::set_value(double val)
{
	if (type != ptDouble)  return false;

	std::ostringstream ss;  ss << val;
	if (ss.fail()) return false;

	_value = ss.str();
	return true;
}

bool TPluginParam::set_value(const std::string & aVal)
{
	if (type != ptString)  return false;

	_value = aVal;
	return true;
}
//-----------------------------------------------------------------------------


bool TPluginParam::get_value(int& val) const
{
	if (type != ptInt)   return false;
	return !(std::istringstream(_value) >> val).fail();
}

bool TPluginParam::get_value(double& val) const
{
	if (type != ptDouble) return false;
	return !(std::istringstream(_value) >> val).fail();
}

bool TPluginParam::get_value(const std::string * &val) const
{
	if (type != ptString) return false;
	val = &_value;  return true;
}
//-----------------------------------------------------------------------------
