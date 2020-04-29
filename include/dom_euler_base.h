#pragma once

#include "mesh.h"

#include "io-field.h"

#include "flux_euler.h"

// TODO: interface for flow model
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
// Domain is something like Scheme, but in a little more general sense
// (experience needed from several flow cases to produce common interface)

// io primitive variables
struct t_PrimVarsIO {

	double u;
	double v;
	double w;
	double p;
	double t;
	t_PrimVarsIO() {}
	t_PrimVarsIO(const t_ConsVars& cv) {

		t_PrimVars pv = cv.calcPrimVars();

		u = pv.getU();
		v = pv.getV();
		w = pv.getW();
		p = pv.getP();
		t = calcTempByRP(pv.getR(), pv.getP());

	}

	t_ConsVars calcConsVars() const {

		t_PrimVars pv;

		pv.setR(calcRhoByPT(this->p,this->t));
		pv.setUVW(t_Vec3(u,v,w));
		pv.setP(p);

		return pv.calcConsVars();
	}

};

class t_DomEuBase : public t_Mesh {

protected:

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

	virtual void allocateFlowSolution();
	virtual void initializeFlow();
	double loadField(std::string fieldName);
	int getNu() const { return NConsVars; }
	std::vector<std::string> getFuncNamesIO() const;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const;


	virtual void makeTimeStep();
	virtual double calcDt() const;
	virtual void calcReconstData() = 0;
	virtual void calcFaceFlux(int iZone, lint iFace) = 0;

	virtual void checkMinMaxCSV();

	// debug
	void dump_flow();
	void dump_geom();
	void checkFlow();
	virtual ~t_DomEuBase();
};

extern t_DomEuBase* G_pDom;
