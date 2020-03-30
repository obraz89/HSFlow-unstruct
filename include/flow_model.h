#pragma once

#include "matrix_small.h"

// TODO: for perfect gas imply that NPrimVars = NConsVars
//static const int NPrimVars = 5;

static const int NConsVars = 5;

using t_VecConsVars = t_Vec<NConsVars>;

struct t_FlowModelParams {

	double Gamma;

	double Pr;

};

extern t_FlowModelParams G_FlowModelParams;

void initialize_flow_model();

class t_ConsVars;

// vector of primitive flow variables 
// (rho, u, v, w, p)
class t_PrimVars : public t_VecConsVars {
public:
	t_PrimVars() = default;
	t_PrimVars(const t_VecConsVars& v) :t_VecConsVars(v) {}
	t_ConsVars calcConsVars() const;
	t_PrimVars& setByCV(const t_ConsVars& cv);
	double getR() const { return data[0]; }
	double getU() const { return data[1]; }
	double getP() const { return data[4]; }
	//double calcRhoE() const;

};

// vector of conservative flow variables
// (rho, rho*u, rho*v, rho*w, rho*E)
class t_ConsVars : public t_VecConsVars {
public:
	t_ConsVars() :t_VecConsVars() {}
	t_ConsVars(const t_VecConsVars& v) :t_VecConsVars(v) {}
	t_ConsVars& setByPV(const t_PrimVars& pv);
	t_PrimVars calcPrimVars() const;
};

class t_Flux : public t_VecConsVars {
public:
	t_Flux() :t_VecConsVars() {}
	t_Flux(const t_VecConsVars& v) :t_VecConsVars(v) {}
	void calc(const t_ConsVars& csv);
	void calc(const t_PrimVars& pvs);

};

// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_Flux& f);

// non-dimensional speed of sound
double calcSoundSpeed(const t_PrimVars& pvs); 