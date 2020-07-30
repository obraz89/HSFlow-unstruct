#pragma once

#include "lsq_uvwpt.h"

#include "dom_ns_base.h"

#include "flux_euler.h"

class t_DomNSLSQ : public t_DomNSBase, public t_LSQData {
public:
	void allocateFlowSolution();
	void calcFaceFlux(int iZone, lint iFace);

	void calcFaceFluxEuler(int iZone, lint iFace, t_FluxEu& fluxEu);
	void calcFaceFluxVisc(int iZone, lint iFace, t_VecConsVars& fluxVisc);

	// implement t_Dom5LSQ
	t_ConsVars calcVirtCellCSV(int iZone, lint iFace) const;

	void calcDataForFaceGradRUVWT(int iZone, lint iFace, t_VecConsVars& Umy, 
		t_VecConsVars& Uop, t_Mat<NConsVars, MaxNumVertsInFace>& UVerts) const;
};
