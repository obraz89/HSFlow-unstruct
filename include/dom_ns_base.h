#pragma once

#include "dom_base_uvwpt.h"

#include "flux_ns.h"

// Flow solution for a zone:
// Fluxes - fluxes at faces
// ConsVars - conservative flow variables at cell centers
// FaceGrdUVWPT - gradients of UVWPT variables at the faces

// IMPORTANT TODO: better storing PVs!!!
struct t_ZoneFlowDataNS {
	t_VecConsVars* Fluxes;
	t_ConsVars* PVCells;
	t_ConsVars* PVVerts;
};

// Domain for Euler equations 
// Domain is something like Scheme, but in a little more general sense
// (experience needed from several flow cases to produce common interface)

class t_DomNSBase : public t_Dom5 {

protected:
	t_ZoneFlowDataNS* ZonesSol;

public:

	t_ConsVars& getCellCSV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].PVCells[cell_id];
	};

	const t_ConsVars& getCellCSV(int zone_id, lint cell_id) const {
		return ZonesSol[zone_id].PVCells[cell_id];
	};

	virtual void allocateFlowSolution();
	virtual void exchangeCSV();

	t_VecConsVars& getFlux(int zone_id, lint face_id) {
		return ZonesSol[zone_id].Fluxes[face_id];
	}

	virtual ~t_DomNSBase();
};

extern t_DomNSBase* G_pDomNS;
