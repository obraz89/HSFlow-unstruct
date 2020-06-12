#include "rs_euler.h"

#include "flow_model_perfect_gas.h"

#include "settings.h"

#include "logging.h"

#include <cmath>

// Davis estimate for wave speeds
// Toro book, Ch.10.5 (p.302)
void calcWaveSpeedDavisEstim(
	const t_PrimVars& pvl, const t_PrimVars& pvr, double& sl, double& sr) {

	t_PrimVars pvm = 0.5 * (pvl + pvr);

	double ul, um, ur;

	ul = pvl.getU();
	um = pvm.getU();
	ur = pvr.getU();

	double cl, cm, cr;

	cl = calcSoundSpeedByRP(pvl.getR(), pvl.getP());
	cm = calcSoundSpeedByRP(pvm.getR(), pvm.getP());
	cr = calcSoundSpeedByRP(pvr.getR(), pvr.getP());

	// sl is minimum of (ul-cl, ur-cr, um-cm)
	sl = fmin(ul - cl, ur - cr);
	sl = fmin(sl, um - cm);

	// sr is maximun of (ul+cl, ur+cr, um+cm)
	sr = fmax(ul + cl, ur + cr);
	sr = fmax(sr, um + cm);

};

double calcWaveSpeedDavisEstimMaxAbs(const t_PrimVars& pvl, const t_PrimVars& pvr) {
	double sl, sr;
	calcWaveSpeedDavisEstim(pvl, pvr, sl, sr);
	return fmax(fabs(sl), fabs(sr));
}

void calcRusanovFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, t_FluxEu& flux) {

	double sl, sr;

	calcWaveSpeedDavisEstim(pvl, pvr, sl, sr);
	double s = fmax(fabs(sl), fabs(sr));

	t_ConsVars cvl, cvr;

	t_FluxEu fl, fr;

	calcCVFlux(pvl, cvl, fl);
	calcCVFlux(pvr, cvr, fr);

	// flux = 0.5 * (fl + fr - s * (cvr - cvl));
	for (int i = 0; i < NConsVars; i++)
		flux[i] = 0.5 * (fl[i] + fr[i] - s * (cvr[i] - cvl[i]));
	

};

void calcRoeFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, t_FluxEu& flux) {

	t_ConsVars cvl, cvr;

	t_FluxEu fl, fr;

	calcCVFlux(pvl, cvl, fl);
	calcCVFlux(pvr, cvr, fr);

	// calculate Roe averages

	double c_denom = 1.0 / (sqrt(pvl.getR()) + sqrt(pvr.getR()));
	double kl = sqrt(pvl.getR()) * c_denom;
	double kr = sqrt(pvr.getR()) * c_denom;

	t_PrimVars pvm;
	{
		double rm = sqrt(pvl.getR() * pvr.getR());
		t_Vec3 Ul = pvl.getUVW();
		t_Vec3 Ur = pvr.getUVW();
		t_Vec3 Um;

		Um = kl * Ul + kr * Ur;

		double Hl = pvl.calcH();
		double Hr = pvr.calcH();
		double Hm = kl * Hl + kr * Hr;
		pvm.setByRhoUH(rm, Um, Hm);
	}

	const double& u = pvm.getU();
	const double& v = pvm.getV();
	const double& w = pvm.getW();

	const double c = calcSoundSpeedByRP(pvm.getR(), pvm.getP());
	const double c_inv = 1.0 / c;
	const double c_inv2 = c_inv * c_inv;

	const double K = G_FlowModelParams.Gamma - 1.0;

	const double q2 = u * u + v * v + w * w;
	// TODO: H: this depends on flow model
	const double H = pvm.calcH();

	t_SqMat<NConsVars> R, Q, L;

	// fill R
	{

		R[0][0] = 1;
		R[0][1] = 1;
		R[0][2] = 1;
		R[0][3] = 0;
		R[0][4] = 0;

		R[1][0] = u - c;
		R[1][1] = u;
		R[1][2] = u + c;
		R[1][3] = 0;
		R[1][4] = 0;

		R[2][0] = v;
		R[2][1] = v;
		R[2][2] = v;
		R[2][3] = 1;
		R[2][4] = 0;

		R[3][0] = w;
		R[3][1] = w;
		R[3][2] = w;
		R[3][3] = 0;
		R[3][4] = 1;

		R[4][0] = H - u * c;
		R[4][1] = 0.5 * q2;
		R[4][2] = H + u * c;
		R[4][3] = v;
		R[4][4] = w;

	}

	// fill Q
	// TODO: entropy fix
	{
		Q[0][0] = fabs(u - c);
		Q[1][1] = fabs(u);
		Q[2][2] = fabs(u + c);
		Q[3][3] = fabs(u);
		Q[4][4] = fabs(u);
	}

	// fill L
	{

		L[0][0] = 0.25 * K * q2 * c_inv2 + 0.5 * u * c_inv;
		L[0][1] = -0.5 * K * u * c_inv2 - 0.5 * c_inv;
		L[0][2] = -0.5 * K * v * c_inv2;
		L[0][3] = -0.5 * K * w * c_inv2;
		L[0][4] = 0.5 * K * c_inv2;

		L[1][0] = 1.0 - 0.5 * K * q2 * c_inv2;
		L[1][1] = K * u * c_inv2;
		L[1][2] = K * v * c_inv2;
		L[1][3] = K * w * c_inv2;
		L[1][4] = -K * c_inv2;

		L[2][0] = 0.25 * K * q2 * c_inv2 - 0.5 * u * c_inv;
		L[2][1] = -0.5 * K * u * c_inv2 + 0.5 * c_inv;
		L[2][2] = -0.5 * K * v * c_inv2;
		L[2][3] = -0.5 * K * w * c_inv2;
		L[2][4] = 0.5 * K * c_inv2;

		L[3][0] = -v;
		L[3][1] = 0;
		L[3][2] = 1;
		L[3][3] = 0;
		L[3][4] = 0;

		L[4][0] = -w;
		L[4][1] = 0;
		L[4][2] = 0;
		L[4][3] = 1;
		L[4][4] = 0;


	}
	// M = R*Q*L
	// TODO: optimize this to reduce number of intermed matrices (here we have 2)
	t_SqMat<NConsVars> M = (R * Q) * L;

	t_ConsVars dfRoe = M * (cvr - cvl);
	// flux = 0.5 * (fl + fr - RQ(|Lambda|)L * (cvr - cvl));
	for (int i = 0; i < NConsVars; i++)
		flux[i] = 0.5 * (fl[i] + fr[i] - dfRoe[i]);

};

void calcRSFlux(const t_PrimVars& pvl, const t_PrimVars& pvr, t_FluxEu& flux) {

	if (g_genOpts.strRiemannSolver.compare("Rusanov") == 0) {

		calcRusanovFlux(pvl, pvr, flux);
		return;

	}

	if (g_genOpts.strRiemannSolver.compare("Roe") == 0) {
		calcRoeFlux(pvl, pvr, flux);
		return;
	}

	hsLogError("calcRSFlux:unknown option for Riemann Solver = %s", g_genOpts.strRiemannSolver.c_str());

};


