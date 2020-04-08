#pragma once

#include "PluginBase.h"

#include <vector>

#include "IniFile.hpp"

#include "common_data.h"

// base class for information of face-type boundary conditions
class t_BCDataFace : public TPlugin{
protected:
	std::string nameOfFldSection;
public:
	t_BCDataFace() = delete;
	t_BCDataFace(const std::string& sect):nameOfFldSection(sect) {}
	const std::string& getSectName() const { return nameOfFldSection; }
	virtual const std::string& getBCKindName() const = 0;
	void init(std::string& ini_data, const std::string& spec) = 0;
	virtual ~t_BCDataFace() {}

};

/**
 * List of bc sets
 // read bc_list section from config dynamically
 // so use iniFile capabilities directly instead of plugin params
 */
class t_BCList : public TPlugin{

	std::map<std::string, t_BCDataFace*>  _pBCs;

public:
	virtual void addBCsetByName(std::string name, std::string bc_kind_str, std::string& ini_data) = 0;
	virtual std::string getSupportedBCsStr() = 0;
	virtual bool has(std::string sectionname) = 0;
	virtual t_FaceBCID getBCID(std::string sectionname) = 0;

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

	virtual ~t_BCList() {
		std::map<std::string, t_BCDataFace*>::iterator it;

		for (auto elem : _pBCs) delete elem.second;

	}

};

extern t_BCList* G_pBCList;