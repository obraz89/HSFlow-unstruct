#include "flow_vars_uvwpt.h"

#include "flow_params.h"

#include "flow_model_perfect_gas.h"

//******************************************t_VecConsVars

void t_VecConsVars::rotate(const t_SqMat3& R) {

	// first component is scalar - not changed

	// components 2-4 rotated as 3d vector
	t_Vec3 res, v;

	for (int i = 0; i < 3; i++) v[i] = data[i + 1];

	res = R * v;

	for (int i = 0; i < 3; i++) data[i + 1] = res[i];

	// fifth component is scalar - not changed

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

t_ConsVars t_PrimVars::calcConsVars() const {

	t_ConsVars csv;

	csv.setByPV(*this);

	return csv;
}

double t_PrimVars::calcVeloSq() const {

	return getU() * getU() + getV() * getV() + getW() * getW();

}

double t_PrimVars::calcH() const {

	double c = calcSoundSpeedByRP(getR(), getP());
	double q2 = calcVeloSq();

	double H = 0.5 * q2 + c * c / (G_FlowModelParams.Gamma - 1.0);
	return H;
};

void t_PrimVars::setByRhoUH(double rho, const t_Vec3& UVW, double H) {

	setR(rho);

	setUVW(UVW);
	double c1 = (G_FlowModelParams.Gamma - 1.0) / G_FlowModelParams.Gamma;
	double p = rho * c1 * (H - 0.5 * calcVeloSq());

	setP(p);

};

void t_PrimVars::setByRhoUT(double rho, const t_Vec3& UVW, double T) {

	setR(rho);

	setUVW(UVW);

	setP(calcPressureByRT(rho, T));

};

void t_PrimVars::setValAtInf() {

	const t_FlowParamsFreeStream& fp = G_FreeStreamParams;

	if (g_genOpts.nonDimType == t_NonDimType::FreeStreamVelo) {

		setR(1.0);

		setUVW(fp.getUInf());

		setP(1.0 / calcGMaMa());

		return;

	}

	if (g_genOpts.nonDimType == t_NonDimType::FreeStreamSoundSpeed) {

		setR(1.0);

		setUVW(fp.getUInf());

		setP(1.0 / G_FlowModelParams.Gamma);

		return;

	}

	hsLogError("t_PrimVars::setValAtInf: unknown non dim option");

}

void t_PrimVars::calcJac(t_SqMat<NConsVars>& Jac) const {

	const double H = calcH();
	const double& Gamma = G_FlowModelParams.Gamma;
	const double q2 = calcVeloSq();

	const double& u = getU();
	const double& v = getV();
	const double& w = getW();

	Jac[0][0] = 0.0;
	Jac[0][1] = 1.0;
	Jac[0][2] = 0.0;
	Jac[0][3] = 0.0;
	Jac[0][4] = 0.0;

	Jac[1][0] = 0.5 * (Gamma - 1.0) * q2 - u * u;
	Jac[1][1] = (3.0 - Gamma) * u;
	Jac[1][2] = (1.0 - Gamma) * v;
	Jac[1][3] = (1.0 - Gamma) * w;
	Jac[1][4] = Gamma - 1.0;

	Jac[2][0] = -1.0 * u * v;
	Jac[2][1] = v;
	Jac[2][2] = u;
	Jac[2][3] = 0.0;
	Jac[2][4] = 0.0;

	Jac[3][0] = -1.0 * u * w;
	Jac[3][1] = w;
	Jac[3][2] = 0.0;
	Jac[3][3] = u;
	Jac[3][4] = 0.0;

	Jac[4][0] = (0.5 * (Gamma - 1.0) * q2 - H) * u;
	Jac[4][1] = H + (1.0 - Gamma) * u * u;
	Jac[4][2] = (1.0 - Gamma) * u * v;
	Jac[4][3] = (1.0 - Gamma) * u * w;
	Jac[4][4] = Gamma * u;

};

double t_PrimVars::calcT() const {
	return calcTempByRP(getR(), getP());
}

t_VecConsVars t_PrimVars::calcUVWPT() const {

	t_VecConsVars UVWPT;

	UVWPT[0] = getU();
	UVWPT[1] = getV();
	UVWPT[2] = getW();
	UVWPT[3] = getP();
	UVWPT[4] = calcT();

	return UVWPT;

};

t_VecConsVars t_PrimVars::calcRUVWT() const{

	t_VecConsVars RUVWT(*this);
	RUVWT[4] = calcT();

	return RUVWT;
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

t_PrimVars t_ConsVars::calcPrimVars() const {

	t_PrimVars pvs;

	pvs.setByCV(*this);

	return pvs;
}


void t_ConsVars::setValAtInf() {

	t_PrimVars prv;

	const t_FlowParamsFreeStream& fp = G_FreeStreamParams;

	prv.setValAtInf();

	*this = prv.calcConsVars();

}

void t_ConsVars::_inflate(const t_SqMat<3>& R, t_SqMat<NConsVars>& dest) {

	dest.reset();

	// first row
	dest[0][0] = 1.0;
	// rows 2-4 from rotation matrix
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dest[i + 1][j + 1] = R[i][j];

	dest[4][4] = 1.0;

}

void t_ConsVars::inflateRotMat(const t_MatRotN& mat_rot_coefs, t_SqMat<NConsVars>& mat) {

	t_SqMat3 R;

	R.set(mat_rot_coefs);

	_inflate(R, mat);

}

void t_ConsVars::inflateRotMatInv(const t_MatRotN& mat_rot_coefs, t_SqMat<NConsVars>& mat) {

	t_SqMat3 R;

	R.set_inv(mat_rot_coefs);

	_inflate(R, mat);

}