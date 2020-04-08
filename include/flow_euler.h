#pragma once

#include "flow_common.h"

#include "flux_euler.h"

struct t_ZoneFlowData {
	// fluxes through
	t_Flux* Fluxes;
	t_ConsVars* ConsVars;
};

class t_DomainEuler : public t_Domain {

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
