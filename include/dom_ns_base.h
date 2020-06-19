#pragma once

#include "dom_base_uvwpt.h"

#include "flux_ns.h"

// TODO: interface for flow model
#include "flow_model_perfect_gas.h"

// Flow solution for a zone:
// Fluxes - fluxes at faces
// ConsVars - conservative flow variables at cell centers
struct t_ZoneFlowDataNS {
	// fluxes through
	t_FluxNS* Fluxes;
	t_ConsVars* ConsVars;
};

// Domain for Euler equations 
// Domain is something like Scheme, but in a little more general sense
// (experience needed from several flow cases to produce common interface)

class t_DomNSBase : public t_Dom5 {

protected:

	t_ZoneFlowDataNS* ZonesSol;

public:

	t_ConsVars& getCellCSV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].ConsVars[cell_id];
	};

	const t_ConsVars& getCellCSV(int zone_id, lint cell_id) const {
		return ZonesSol[zone_id].ConsVars[cell_id];
	};

	void allocateFlowSolution();
	virtual void exchangeCSV();

	t_FluxNS& getFlux(int zone_id, lint face_id) {
		return ZonesSol[zone_id].Fluxes[face_id];
	}

	t_VecConsVars& getFlux5(int zone_id, lint face_id) {
		return getFlux(zone_id, face_id);
	}

	virtual ~t_DomNSBase();
};

extern t_DomNSBase* G_pDomNS;
