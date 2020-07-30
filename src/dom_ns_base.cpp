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
		lint nVerts = zne.getnVerts();

		fdata.Fluxes = new t_VecConsVars[nFaces];
		fdata.PVCells = new t_ConsVars[nCellsTot];
		fdata.PVVerts = new t_ConsVars[nVerts];


	}

}

void t_DomNSBase::initializeFlow() {

	// read CSV, exchange CSV
	t_Dom5::initializeFlow();

	// calculate vertex values via cell csvs
	calcVertexValues();
}

void t_DomNSBase::makeTimeStep() {

	// basic explicit step, update CSV of cells
	t_Dom5::makeTimeStep();

	// update vertex values
	calcVertexValues();

}

void t_DomNSBase::exchangeCSV() { G_GhostMngNS.exchangeCSV(); }

t_DomNSBase::~t_DomNSBase() {

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowDataNS& fdata = ZonesSol[i];

		delete[] fdata.Fluxes;
		delete[] fdata.PVCells;
		delete[] fdata.PVVerts;

	}

	delete[] ZonesSol;

}

void t_DomNSBase::calcCellWeightsForVertices() {

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];

		for (int ivert = 0; ivert < zne.getnVerts(); ivert++) {

			zne.getVert(ivert).calcAllocNeigCoefs();

		}
	}

};

void t_DomNSBase::calcFaceGradMatrices() {

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];

		for (int iFace = 0; iFace < zne.getNFaces(); iFace++) {

			zne.getFace(iFace).ComputeMatGrad();

		}
	}

}

void t_DomNSBase::calcVertexValues() {


	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_ConsVars csv_vert;
		t_ConsVars csv_neig;

		for (int iVert = 0; iVert < zne.getnVerts(); iVert++) {

			csv_vert.reset();

			t_Vert& vert = zne.getVert(iVert);

			for (int j = 0; j < vert.NNeigCells; j++) {

				csv_neig = getCellCSV(iZone, vert.pNeigCells[j]->Id);

				for (int k = 0; k < NConsVars; k++)
					csv_vert[k] += vert.pNeigCoefs[j] * csv_neig[k];

			}
				
			getVertCSV(iZone, iVert) = csv_vert;

		}
	}

};