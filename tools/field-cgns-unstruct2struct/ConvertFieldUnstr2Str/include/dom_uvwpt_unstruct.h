#pragma once

#include "dom_base_unstruct.h"

#include "io_field_unstruct.h"

static const int NConsVars = 5;

struct t_PrimVarsIO {

	double u;
	double v;
	double w;
	double p;
	double t;

};

struct t_ZoneFlowData {
	t_PrimVarsIO* PV = nullptr;
};

class t_Dom5 : public t_DomBase {
protected:
	t_ZoneFlowData* ZonesSol = nullptr;

public:

	t_PrimVarsIO& getCellPV(int zone_id, lint cell_id) {
		return ZonesSol[zone_id].PV[cell_id];
	};

	const t_PrimVarsIO& getCellPV(int zone_id, lint cell_id) const {
		return ZonesSol[zone_id].PV[cell_id];
	};

	virtual void allocateFlowSolution();
	virtual void initializeFlow();
	double loadField(std::string fieldName);
	int getNu() const { return NConsVars; }
	std::vector<std::string> getFuncNamesIO() const;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const;

	virtual ~t_Dom5();

	// dummy implementation of domain base
	virtual void prepareBeforeTimeMarch() {};
	virtual void checkMinMaxCSV() {};
	virtual void makeTimeStep() {};
};


