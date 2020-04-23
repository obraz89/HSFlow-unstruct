#include "dom_euler_lsq.h"

void t_DomEuLSQ::allocateFlowSolution() {

	t_DomEuBase::allocateFlowSolution();

	ZonesRecData = new t_ZoneReconstData[nZones];

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneReconstData& fdata = ZonesRecData[i];

		lint nCellsTot = zne.getnCellsTot();

		fdata.ReconstData = new t_ReconstDataLSQ[nCellsTot];


	}
}

t_DomEuLSQ::~t_DomEuLSQ() {

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneReconstData& fdata = ZonesRecData[i];

		delete[] fdata.ReconstData;

	}

	delete[] ZonesRecData;

}

void t_DomEuLSQ::calcFaceFlux(int iZone, lint iFace) {

	hsLogError("Not Implemented");

};