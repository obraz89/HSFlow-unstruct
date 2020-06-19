#include "flux_ns.h"

#include "flow_model_perfect_gas.h"


//******************************************Flux
// we are in some reference frame (x,y,z) which is rotation from global (X,Y,Z)
// compute inviscid face flux through area with normal (1;0;0) 
// velocities must be in reference frame (x,y,z)!
void t_FluxNS::calc(const t_PrimVars& pv) {

	const double& r = pv[0];
	const double& u = pv[1];
	const double& v = pv[2];
	const double& w = pv[3];
	const double& p = pv[4];

	double rhoE = 0.5 * r * (u * u + v * v + w * w) + p / (G_FlowModelParams.Gamma - 1);

	data[0] = r * u;
	data[1] = r * u * u + p;
	data[2] = r * u * v;
	data[3] = r * u * w;
	data[4] = u * (rhoE + p);

}

// TODO: avoid this function
// use faster calcCVFlux()
void t_FluxNS::calc(const t_ConsVars& cv) {

	calc(cv.calcPrimVars());

}

// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_FluxNS& f) {

	cv.setByPV(pv);
	f.calc(pv);

};
