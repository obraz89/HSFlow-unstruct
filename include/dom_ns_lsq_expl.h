#pragma once

// get reconstr types
#include "dom_euler_lsq.h"

#include "dom_ns_base.h"

class t_DomNSLSQ : public t_DomNSBase, public t_LSQData {
public:
	void allocateFlowSolution();
	void calcFaceFlux(int iZone, lint iFace);

	// implement t_Dom5LSQ
	t_ConsVars calcVirtCellCSV(int iZone, lint iFace) const;
};
