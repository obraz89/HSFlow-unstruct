#pragma once

#include "PluginBase.h"

#include <vector>

#include "IniFile.hpp"

#include "common_data.h"

// base class for information of face-type boundary conditions
class t_BCDataFace : public TPlugin{
protected:
	std::string FamilyName;
public:
	t_BCDataFace() = delete;
	t_BCDataFace(const std::string& a_name):FamilyName(a_name) {}
	const std::string& getFamName() const { return FamilyName; }
	virtual std::string getBCKindName() const = 0;
	// implement TPlugin
	std::string get_name() const { return getBCKindName() + "/" + FamilyName; }
	virtual void init(std::string& ini_data, const std::string& spec) { TPlugin::init(ini_data, spec); };
	virtual void default_settings() {}
	virtual ~t_BCDataFace() {}

};

/**
 * List of bc sets
 // read bc_list section from config dynamically
 // so use iniFile capabilities directly instead of plugin params
 */
class t_BCList : public TPlugin{

public:
	virtual void addBCsetByName(std::string name, std::string bc_kind_str, std::string& ini_data) = 0;
	virtual std::string getSupportedBCsStr() const = 0;
	virtual bool has(std::string sectionname) const = 0;
	virtual t_FaceBCID getBCID(std::string fam_name) const = 0;

	virtual void init(std::string& ini_data, const std::string& spec) {

		TPlugin::init(ini_data, spec);

		ini::IniFile ini;   ini.decode(ini_data);
		std::string bc_list_sname = get_name();
		if (!ini.has(bc_list_sname)) {
			ini::IniSection sect;
			// add some defaults
			sect["Supported_bc_names"] = getSupportedBCsStr();
			sect["inflow"] = "bc_inflow_of_your_model";
			sect["out"] = "bc_outflow_of_your_model";

			ini[bc_list_sname] = sect;
		}

		ini_data = ini.encode();

		for (const auto& field_pair : ini[bc_list_sname]) {

			if (field_pair.first.compare("Supported_bc_names") == 0) {
				// skip this helper option
				continue;
			}

			std::string bc_set_name = field_pair.first;
			std::string bc_kind_name = field_pair.second.asString();

			addBCsetByName(bc_set_name, bc_kind_name, ini_data);

		}
	};

	virtual ~t_BCList() {}

};

extern t_BCList* G_pBCList;