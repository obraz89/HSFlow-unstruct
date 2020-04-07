#include "flow_model.h"

#include "flow_params.h"

t_FlowModelParams G_FlowModelParams;

void initialize_flow_model() {

	G_FlowModelParams.Gamma = 1.4;

	G_FlowModelParams.Pr = 0.72;

}

//******************************************t_VecConsVars

t_VecConsVars& t_VecConsVars::rotate(const t_SqMat3& R) {

	// first component is scalar - not changed

	// components 2-4 rotated as 3d vector
	t_Vec3 res, v;

	for (int i = 0; i < 3; i++) v[i] = data[i+1];

	res = R * v;

	for (int i = 0; i < 3; i++) data[i+1] = res[i];

	// fifth component is scalar - not changed

	return *this;

}

//******************************************PrimVars

t_PrimVars& t_PrimVars::setByCV(const t_ConsVars& cv) {

	const double& r = cv[0];

	double r_inv = 1.0 / r;

	double u = r_inv * cv[1];
	double v = r_inv * cv[2];
	double w = r_inv * cv[3];

	double p = (G_FlowModelParams.Gamma - 1) * (cv[4] - 0.5 * r * (u * u + v * v + w * w));

	data[0] = r;
	data[1] = u;
	data[2] = v;
	data[3] = w;
	data[4] = p;

	return *this;
}

t_ConsVars t_PrimVars::calcConsVars() const{

	t_ConsVars csv;
	
	csv.setByPV(*this);

	return csv;
}
//******************************************ConsVars

t_ConsVars& t_ConsVars::setByPV(const t_PrimVars& pv) {

	const double& r = pv[0];
	const double& u = pv[1];
	const double& v = pv[2];
	const double& w = pv[3];
	const double& p = pv[4];

	data[0] = r;

	data[1] = r * u;

	data[2] = r * v;

	data[3] = r * w;

	data[4] = 0.5 * r * (u * u + v * v + w * w) + p / (G_FlowModelParams.Gamma - 1);

	return *this;

}

t_PrimVars t_ConsVars::calcPrimVars() const{

	t_PrimVars pvs;

	pvs.setByCV(*this);

	return pvs;
}

//******************************************Flux
// we are in some reference frame (x,y,z) which is rotation from global (X,Y,Z)
// compute inviscid face flux through area with normal (1;0;0) 
// velocities must be in reference frame (x,y,z)!
void t_Flux::calc(const t_PrimVars& pv) {

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
void t_Flux::calc(const t_ConsVars& cv) {

	calc(cv.calcPrimVars());

}

// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_Flux& f) {

	cv.setByPV(pv);
	f.calc(pv);

};

// common functions

double calcSoundSpeed(const t_PrimVars& pvs) {

	return sqrt(G_FlowModelParams.Gamma * pvs.getP() / pvs.getR());

}

double calcGMaMaInv() {
	double M2 = G_FreeStreamParams.getMach() * G_FreeStreamParams.getMach();
	return 1.0 / (G_FlowModelParams.Gamma * M2);
}

t_ConsVars calcConsVarsInf() {

	t_PrimVars prv;

	const t_FlowParamsFreeStream& fp = G_FreeStreamParams;

	prv.setR(1.0);

	prv.setUVW(fp.getUInf());

	prv.setP(calcGMaMaInv());

	t_ConsVars csv = prv.calcConsVars();

	return csv;

}
