#pragma once

#include "dom_euler_base.h"
#include "lsq_uvwpt.h"

class t_DomEuLSQ : public t_DomEuBase, public t_LSQData {

public:
	void allocateFlowSolution();
	void calcFaceFlux(int iZone, lint iFace);

	//implement t_DomBase
	void prepareBeforeTimeMarch();

	// implement t_Dom5LSQ
	t_ConsVars calcVirtCellCSV(int iZone, lint iFace) const;
};

