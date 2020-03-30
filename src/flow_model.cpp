#include "flow_model.h"

t_FlowModelParams G_FlowModelParams;

void initialize_flow_model() {

	G_FlowModelParams.Gamma = 1.4;

	G_FlowModelParams.Pr = 0.72;

}

//******************************************PrimVars

t_ConsVars t_PrimVars::calcConsVars() const{

	t_ConsVars csv;

	const double& r = data[0];
	const double& u = data[1];
	const double& v = data[2];
	const double& w = data[3];
	const double& p = data[4];

	csv[0] = r;
	
	csv[1] = r * u;
	
	csv[2] = r * v;
	
	csv[3] = r * w;
	
	csv[4] = 0.5*r*(u*u + v*v + w*w) + p / (G_FlowModelParams.Gamma - 1);

	return csv;
}
//******************************************ConsVars

t_PrimVars t_ConsVars::calcPrimVars() const{

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

//******************************************Flux
// we are in some reference frame (x,y,z) which is rotation from global (X,Y,Z)
// compute inviscid face flux through area with normal (1;0;0) 
// velocities must be in reference frame (x,y,z)!
void t_Flux::computeFlux(const t_PrimVars& pvs) {

	const double& r = data[0];
	const double& u = data[1];
	const double& v = data[2];
	const double& w = data[3];
	const double& p = data[4];

	double rhoE = 0.5 * r * (u * u + v * v + w * w) + p / (G_FlowModelParams.Gamma - 1);

	data[0] = r * u;
	data[1] = r * u * u + p;
	data[2] = r * u * v;
	data[3] = r * u * w;
	data[4] = u * (rhoE + p);

}