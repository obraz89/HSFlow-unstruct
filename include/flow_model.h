#pragma once

//#include "common_data.h"

// TODO: for perfect gas imply that NPrimVars = NConsVars
//static const int NPrimVars = 5;

static const int NConsVars = 5;

struct t_FlowModelParams {

	double Gamma;

	double Pr;

};

extern t_FlowModelParams G_FlowModelParams;

void initialize_flow_model();

// vector of size NConsVars
class t_VecConsVars {
protected:
	double data[NConsVars];
public:
	double& operator[](int ind) { return data[ind]; }
	const double& operator[](int ind) const { return data[ind]; }

};

class t_ConsVars;

// vector of primitive flow variables 
// (rho, u, v, w, p)
class t_PrimVars : public t_VecConsVars {
public:
	t_ConsVars toConsVars();

};

// vector of conservative flow variables
// (rho, rho*u, rho*v, rho*w, rho*E)
class t_ConsVars : public t_VecConsVars {
public:
	t_PrimVars toPrimVars();
};
