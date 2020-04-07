#include "common_data.h"

#include "flow_model.h"

#include "ghost_manager.h"

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