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

// virtual cell abutts real cell via BC face
struct t_ZoneVirtCells {
	t_Cell* VirtCells;
};

class t_DomEuLSQ : public t_DomEuBase {
	t_ZoneReconstData* ZonesRecData;
	t_ZoneVirtCells* ZonesVirtCells;
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
	void calcCellGradCSV(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGrad);
	void calcReconstDataLSQ(int iZone, lint iCell);

	void initializeVirtCells();

	t_ConsVars calcVirtCellCSV(int iZone, lint iFace) const;
	t_ConsVars getNeigCellCSV(int iZone, lint iCell, int indFace) const;

	t_Vec<NConsVars> calcSlopeLimiters(int iZone, lint iCell, const t_Mat<NConsVars, 3>& CellGradPV) const;
	~t_DomEuLSQ();
};

