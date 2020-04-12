#pragma once

#include "mesh.h"

#include "flux_euler.h"

#include "io-field.h"

#include "flow_model_perfect_gas.h"

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

// io primitive variables

struct t_PrimVarsIO {

	double u;
	double v; 
	double w;
	double p; 
	double t;

	t_PrimVarsIO(const t_ConsVars& cv) {

		t_PrimVars pv = cv.calcPrimVars();

		u = pv.getU();
		v = pv.getV();
		w = pv.getW();
		p = pv.getP();
		t = calcTempByRP(pv.getR(), pv.getP());

	}

};

class t_DomainEuler : public t_Mesh {

	t_ZoneFlowData* ZonesSol;

public:

	t_ConsVars& getCellCSV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].ConsVars[cell_id];
	};

	const t_ConsVars& getCellCSV(int zone_id, lint cell_id) const{
		return ZonesSol[zone_id].ConsVars[cell_id];
	};

	t_Flux& getFlux(int zone_id, lint face_id) {
		return ZonesSol[zone_id].Fluxes[face_id];
	}

	// implement domain virtuals
	void allocateFlowSolution();
	void initializeFlow();
	double loadField(std::string fieldName);
	int getNu() const { return NConsVars; }
	std::vector<std::string> getFuncNamesIO() const;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const;


	void makeTimeStep();
	double calcDt();
	void calcFaceFlux(int iZone, lint iFace);
	void dump_flow();
	void dump_geom();
	~t_DomainEuler();
};

extern t_DomainEuler G_Domain;
