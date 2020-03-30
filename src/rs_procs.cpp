#include "rs_procs.h"

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

	cl = calcSoundSpeed(pvl);
	cm = calcSoundSpeed(pvm);
	cr = calcSoundSpeed(pvr);

	// sl is minimum of (ul-cl, ur-cr, um-cm)
	sl = fmin(ul - cl, ur - cr);
	sl = fmin(sl, um - cm);

	// sr is maximun of (ul+cl, ur+cr, um+cm)
	sr = fmax(ul + cl, ur + cr);
	sr = fmax(sr, um + cm);

};

void calcRusanovFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, t_Flux& flux) {

	double sl, sr;

	calcWaveSpeedDavisEstim(pvl, pvr, sl, sr);
	double s = fmax(fabs(sl), fabs(sr));

	t_ConsVars cvl, cvr;

	t_Flux fl, fr;

	calcCVFlux(pvl, cvl, fl);
	calcCVFlux(pvr, cvr, fr);

	// flux = 0.5 * (fl + fr - s * (cvr - cvl));
	for (int i = 0; i < NConsVars; i++)
		flux[i] = 0.5 * (fl[i] + fr[i] - s * (cvr[i] - cvl[i]));
	

};


