#include "flow_euler.h"

#include "ghost_manager.h"

#include <fstream>


t_DomainEuler G_Domain;

void t_DomainEuler::allocateFlowSolution() {

	ZonesSol = new t_ZoneFlowData[nZones];

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowData& fdata = ZonesSol[i];
		
		lint nFaces = zne.getNFaces();
		lint nCellsTot = zne.getnCellsTot();

		fdata.Fluxes = new t_Flux[nFaces];
		fdata.ConsVars = new t_ConsVars[nCellsTot];


	}

}

void t_DomainEuler::initializeFlow() {

	// set real cell values
	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		for (int i = 0; i < zne.getnCellsReal(); i++) {

			pcell = zne.getpCell(i);
			getCellCSV(iZone, i).setValAtInf();
		}

	}
	// set ghost values
	G_GhostManager.exchangeCSV();

}

void t_DomainEuler::dump_flow() {

	std::string fn("dump_flow.txt");

	std::ofstream ofstr(fn);

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		ofstr << "=========Zone #" << iZone << "===========\n";

		for (int i = 0; i < zne.getnCellsTot(); i++) {

			pcell = zne.getpCell(i);
			ofstr << "cell #" << i << getCellCSV(iZone, i).to_str();
		}

	}

	ofstr.flush();

}

void t_DomainEuler::dump_geom() {



}

t_DomainEuler::~t_DomainEuler() {

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowData& fdata = ZonesSol[i];

		lint nFaces = zne.getNFaces();
		lint nCellsTot = zne.getnCellsTot();

		delete[] fdata.Fluxes;
		delete[] fdata.ConsVars;

	}

	delete[] ZonesSol;

}