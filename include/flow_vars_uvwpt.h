#pragma once

#include "matrix_small.h"

static const int NConsVars = 5;

class t_VecConsVars : public t_Vec<NConsVars> {

public:
	t_VecConsVars() :t_Vec<NConsVars>() {}
	t_VecConsVars(const t_Vec<NConsVars>& v) : t_Vec<NConsVars>(v) {}
	// Conservative & Primitive vectors rotate the same way
	virtual void rotate(const t_SqMat3& R);
	virtual std::string to_str() const { return t_Vec<NConsVars>::to_str(); }

};

class t_ConsVars;

// vector of primitive flow variables 
// (rho, u, v, w, p)
// for perfect gas imply that NPrimVars = NConsVars
// for this reason PrimVars inherit from t_VecConsVars
class t_PrimVars : public t_VecConsVars {
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

	double calcH() const;
	double calcVeloSq() const;

	t_Vec3 getUVW() const { return t_Vec3(getU(), getV(), getW()); }

	double getP() const { return data[4]; }
	void setP(double val) { data[4] = val; }

	void setByRhoUH(double rho, const t_Vec3& UVW, double H);
	void setByRhoUT(double rho, const t_Vec3& UVW, double T);

	void setValAtInf();

	void calcJac(t_SqMat<NConsVars>& Jac) const;

	double calcT() const;

	t_VecConsVars calcUVWPT() const;
	t_VecConsVars calcRUVWT() const;

};

// vector of conservative flow variables
// (rho, rho*u, rho*v, rho*w, rho*E)
class t_ConsVars : public t_VecConsVars {
	static void _inflate(const t_SqMat<3>& R, t_SqMat<NConsVars>& dest);
public:
	t_ConsVars() :t_VecConsVars() {}
	t_ConsVars(const t_VecConsVars& v) :t_VecConsVars(v) {}
	t_ConsVars(const t_Vec<NConsVars>& v) : t_VecConsVars(v) {}
	t_ConsVars& setByPV(const t_PrimVars& pv);
	t_PrimVars calcPrimVars() const;

	void setValAtInf();

	static void inflateRotMat(const t_MatRotN& mat_rot, t_SqMat<NConsVars>& mat);
	static void inflateRotMatInv(const t_MatRotN& mat_rot, t_SqMat<NConsVars>& mat);
};


