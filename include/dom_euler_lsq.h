#pragma once

#include "dom_euler_base.h"

struct t_ReconstDataLSQ {

	// scale used to make reconstruction matrix O(1)
	double r;
	// inverse reconst matrix multiplied by r*r
	t_SqMat3 MInvRR;

};

struct t_ZoneReconstData {
	t_ReconstDataLSQ* ReconstData;
};

class t_DomEuLSQ : public t_DomEuBase {
	t_ZoneReconstData* ZonesRecData;
public:
	void allocateFlowSolution();
	void calcReconstData();
	void calcFaceFlux(int iZone, lint iFace);
	// internals
	t_ReconstDataLSQ& getReconstData(int iZone, lint iCell) {
		return ZonesRecData[iZone].ReconstData[iCell];
	};
	const t_ReconstDataLSQ& getReconstData(int iZone, lint iCell) const{
		return ZonesRecData[iZone].ReconstData[iCell];
	};
	void calcCellGradPrimVars(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGrad);
	void calcReconstDataLSQ(int iZone, lint iCell);
	~t_DomEuLSQ();
};

