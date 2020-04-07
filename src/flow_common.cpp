#include "common_data.h"

#include "flow_model.h"

#include "ghost_manager.h"

#include <fstream>

void t_Domain::initializeFlow() {

	t_ConsVars cvs = calcConsVarsInf();

	// set real cell values
	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		for (int i = 0; i < zne.getnCellsReal(); i++) {

			pcell = zne.getpCell(i);
			pcell->ConsVars = cvs;
		}

	}
	// set ghost values
	G_GhostManager.exchangeCSV();

}

void t_Domain::dump_flow() {

	std::string fn("dump_flow.txt");

	std::ofstream ofstr(fn);

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		ofstr << "=========Zone #" << iZone << "===========\n";

		for (int i = 0; i < zne.getnCellsTot(); i++) {

			pcell = zne.getpCell(i);
			ofstr <<"cell #"<<i << pcell->ConsVars.to_str();
		}

	}

	ofstr.flush();

}

void t_Domain::dump_geom() {



}