#include "settings.h"


void load_settings(const std::string& fn) {

	ini::IniFile ini(fn);

	G_Domain.nDim = 3;

	for (auto& sect : ini) {

		hsLogMessage("Parsing ini: has section %s", &sect.first[0]);

	}

};