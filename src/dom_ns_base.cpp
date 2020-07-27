#include "mpi.h"

#include "dom_ns_base.h"

#include "ghost_ns.h"

#include "flux_ns.h"

// TODO: make interface for flow model
#include "flow_model_perfect_gas.h"

#include "settings.h"

#include "common_data.h"

#include <fstream>

void t_DomNSBase::allocateFlowSolution() {

	ZonesSol = new t_ZoneFlowDataNS[nZones];

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowDataNS& fdata = ZonesSol[i];

		lint nFaces = zne.getNFaces();
		lint nCellsTot = zne.getnCellsTot();

		fdata.Fluxes = new t_VecConsVars[nFaces];
		fdata.ConsVars = new t_ConsVars[nCellsTot];
		fdata.FaceGrdUVWPT = new t_Mat<NConsVars, 3>[nFaces];


	}

}

void t_DomNSBase::exchangeCSV() { G_GhostMngNS.exchangeCSV(); }

t_DomNSBase::~t_DomNSBase() {

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowDataNS& fdata = ZonesSol[i];

		delete[] fdata.Fluxes;
		delete[] fdata.ConsVars;

	}

	delete[] ZonesSol;

}