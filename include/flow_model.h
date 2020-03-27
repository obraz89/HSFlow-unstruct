#pragma once

//#include "common_data.h"

static const int NPrimVars = 5;

static const int NConsVars = 5;

struct t_FlowModelParams {

	double Gamma;

	double Pr;

};

extern t_FlowModelParams G_FlowModelParams;

void initialize_flow_model();

class t_ConsVars;

// vector of primitive flow variables 
// (rho, u, v, w, p)
class t_PrimVars {
	double data[NPrimVars];
public:
	t_ConsVars toConsVars();
	double& operator[](int ind) { return data[ind]; }
	const double& operator[](int ind) const{ return data[ind]; }

};

// vector of conservative flow variables
// (rho, rho*u, rho*v, rho*w, rho*E)
class t_ConsVars {
	double data[NConsVars];
public:
	t_PrimVars toPrimVars();
	double& operator[](int ind) { return data[ind]; }
	const double& operator[](int ind) const { return data[ind]; }
};
