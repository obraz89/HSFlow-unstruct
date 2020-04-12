#pragma once

#include "matrix_small.h"

static const int NConsVars = 5;

class t_VecConsVars : public t_Vec<NConsVars> {

public:
	t_VecConsVars() :t_Vec<NConsVars>() {}
	t_VecConsVars(const t_Vec<NConsVars>& v) : t_Vec<NConsVars>(v) {}
	// Conservative & Primitive vectors rotate the same way
	virtual void rotate(const t_SqMat3& R);
	virtual std::string to_str() const {return t_Vec<NConsVars>::to_str(); }

};

class t_ConsVars;

// vector of primitive flow variables 
// (rho, u, v, w, p)
// for perfect gas imply that NPrimVars = NConsVars
// for this reason PrimVars inherit from t_VecConsVars
class t_PrimVars : public t_VecConsVars  {
public:
	t_PrimVars() :t_VecConsVars() {};
	t_PrimVars(const t_VecConsVars& v) :t_VecConsVars(v) {}
	t_PrimVars(const t_Vec<NConsVars>& v) :t_VecConsVars(v) {}
	t_ConsVars calcConsVars() const;
	t_PrimVars& setByCV(const t_ConsVars& cv);

	double getR() const { return data[0]; }
	void setR(double val) { data[0] = val; }

	double getU() const { return data[1]; }
	void setU(double val) { data[1] = val; }

	double getV() const { return data[2]; }
	double getW() const { return data[3]; }

	void setUVW(const t_Vec3& v) {
		for (int i = 0; i < 3; i++) data[i + 1] = v[i];
	}
	double getP() const { return data[4]; }
	void setP(double val) { data[4] = val; }
	//double calcRhoE() const;

	void setValAtInf();

};

// vector of conservative flow variables
// (rho, rho*u, rho*v, rho*w, rho*E)
class t_ConsVars : public t_VecConsVars {
public:
	t_ConsVars() :t_VecConsVars() {}
	t_ConsVars(const t_VecConsVars& v) :t_VecConsVars(v) {}
	t_ConsVars& setByPV(const t_PrimVars& pv);
	t_PrimVars calcPrimVars() const;

	void setValAtInf();
};

class t_Flux : public t_VecConsVars{
public:
	t_Flux() :t_VecConsVars() {}
	t_Flux(const t_VecConsVars& v) :t_VecConsVars(v) {}
	void calc(const t_ConsVars& csv);
	void calc(const t_PrimVars& pvs);

};

// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_Flux& f);


