#pragma once

#include "dom_base.h"

#include "io-field.h"

#include "flow_vars_uvwpt.h"

#include "flow_model_perfect_gas.h"

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

		pv.setR(calcRhoByPT(this->p, this->t));
		pv.setUVW(t_Vec3(u, v, w));
		pv.setP(p);

		return pv.calcConsVars();
	}

};

class t_Dom5 : public t_DomBase {

public:
	virtual t_ConsVars& getCellCSV(int zone_id, lint cell_id) = 0;
	virtual const t_ConsVars& getCellCSV(int zone_id, lint cell_id) const = 0;

	virtual void allocateFlowSolution() = 0;
	virtual void initializeFlow();
	double loadField(std::string fieldName);
	int getNu() const { return NConsVars; }
	std::vector<std::string> getFuncNamesIO() const;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const;

	virtual double calcDt() const;
	virtual void calcFaceFlux(int iZone, lint iFace) = 0;

	virtual void checkMinMaxCSV();
	virtual void makeTimeStep();

	virtual void exchangeCSV() = 0;
	virtual t_VecConsVars& getFlux5(int zone_id, lint face_id) = 0;
};


