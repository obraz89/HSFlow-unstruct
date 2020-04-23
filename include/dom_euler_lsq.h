#pragma once

#include "dom_euler_base.h"

#include "reconstruct_lsq.h"

struct t_ZoneReconstData {
	t_ReconstDataLSQ* ReconstData;
};

class t_DomEuLSQ : public t_DomEuBase {
	t_ZoneReconstData* ZonesRecData;
public:
	void allocateFlowSolution();
	void calcFaceFlux(int iZone, lint iFace);
	~t_DomEuLSQ();
};

