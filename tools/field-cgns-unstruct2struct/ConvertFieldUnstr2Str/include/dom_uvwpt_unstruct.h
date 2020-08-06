#pragma once

#include "dom_base_unstruct.h"

#include "io_field_unstruct.h"

struct t_PrimVars : public t_Vec<NConsVars> {

	double getU() const{ return data[0]; }
	double& getU() { return data[0]; }

	double getV() const{ return data[1]; }
	double& getV() { return data[1]; }

	double getW() const{ return data[2]; }
	double& getW() { return data[2]; }

	double getP() const{ return data[3]; }
	double& getP() { return data[3]; }

	double getT() const{ return data[4]; }
	double& getT() { return data[4]; }

};

struct t_ZoneFlowData {
	t_PrimVars* PV = nullptr;
};

class t_Dom5 : public t_DomBase {
protected:
	t_ZoneFlowData* ZonesSol = nullptr;

public:

	t_PrimVars& getCellPV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].PV[cell_id];
	};

	const t_PrimVars& getCellPV(int zone_id, lint cell_id) const {
		return ZonesSol[zone_id].PV[cell_id];
	};

	virtual void allocateFlowSolution();
	virtual void initializeFlow();
	double loadField(std::string fieldName);
	int getNu() const { return NConsVars; }
	std::vector<std::string> getFuncNamesIO() const;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const;

	t_PrimVars calcVertexPV(int iZone, int iVert) const;

	virtual ~t_Dom5();

	// dummy implementation of domain base
	virtual void prepareBeforeTimeMarch() {};
	virtual void checkMinMaxCSV() {};
	virtual void makeTimeStep() {};
};

extern t_Dom5 G_DomUnst;


