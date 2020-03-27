#include "flow_model.h"

t_FlowModelParams G_FlowModelParams;

void initialize_flow_model() {

	G_FlowModelParams.Gamma = 1.4;

	G_FlowModelParams.Pr = 0.72;

}

t_ConsVars t_PrimVars::toConsVars() {

	t_ConsVars csv;

	double r = data[0];
	double u = data[1];
	double v = data[2];
	double w = data[3];
	double p = data[4];

	csv[0] = r;
	
	csv[1] = r * u;
	
	csv[2] = r * v;
	
	csv[3] = r * w;
	
	csv[4] = 0.5*r*(u*u + v*v + w*w) + p / (G_FlowModelParams.Gamma - 1);

	return csv;
}

t_PrimVars t_ConsVars::toPrimVars() {

	t_PrimVars pvs;

	double r = data[0];
	double u = data[1] / r;
	double v = data[2] / r;
	double w = data[3] / r;
	double p = (G_FlowModelParams.Gamma - 1)*(data[4] - 0.5*r*(u*u + v*v + w*w));

	pvs[0] = r;
	pvs[1] = u;
	pvs[2] = v;
	pvs[3] = w;
	pvs[4] = p;

	return pvs;
}