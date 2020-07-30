#pragma once

#pragma once

#include "dom_base_uvwpt.h"

#include "ghost_uvwpt.h"

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

// shared lsq methods
class t_LSQData {
	t_Dom5* Dom;
	t_ZoneReconstData* ZonesRecData;
	t_ZoneVirtCells* ZonesVirtCells;
protected:
	void allocateLSQData(t_Dom5* a_dom);
	void calcReconstData(t_GhostMng5& ghost_mng);

	// internals
	t_ReconstDataLSQ& getReconstData(int iZone, lint iCell) {
		return ZonesRecData[iZone].ReconstData[iCell];
	};
	const t_ReconstDataLSQ& getReconstData(int iZone, lint iCell) const {
		return ZonesRecData[iZone].ReconstData[iCell];
	};
	void calcCellGradCSV(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGrad);

	// mainly for debug & comparison with face grad data
	void calcCellGradRUVWT(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGrad);
	void calcReconstDataLSQ(int iZone, lint iCell);

	t_Vec<NConsVars> calcSlopeLimiters(int iZone, lint iCell, const t_Mat<NConsVars, 3>& CellGradPV) const;

	void initializeVirtCells();

	virtual t_ConsVars getNeigCellCSV(int iZone, lint iCell, int indFace) const;

	// pure virtuals
	virtual t_ConsVars calcVirtCellCSV(int iZone, lint iFace) const = 0;

	~t_LSQData();
};


