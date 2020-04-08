#pragma once

#include "mesh.h"

#include "flux_euler.h"

// Flow solution for a zone:
// Fluxes - fluxes at faces
// ConsVars - conservative flow variables at cell centers
struct t_ZoneFlowData {
	// fluxes through
	t_Flux* Fluxes;
	t_ConsVars* ConsVars;
};

// Domain for Euler equations
// TODO: inherit from common t_Domain
// (experience needed from several flow cases to produce common interface)

class t_DomainEuler : public t_Mesh {

	t_ZoneFlowData* ZonesSol;

public:

	t_ConsVars& getCellCSV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].ConsVars[cell_id];
	};

	t_Flux& getFlux(int zone_id, lint face_id) {
		return ZonesSol[zone_id].Fluxes[face_id];
	}

	void allocateFlowSolution();
	void initializeFlow();
	void makeTimeStep();
	double calcDt();
	void calcFaceFlux(int iZone, lint iFace);
	void dump_flow();
	void dump_geom();
	~t_DomainEuler();
};

extern t_DomainEuler G_Domain;
