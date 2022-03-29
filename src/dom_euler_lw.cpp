#include "dom_euler_lw.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

void calcLaxWendroffFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, const double dt_dx, t_FluxEu& flux);

void t_DomEuLW::calcFaceFlux(int iZone, lint iFace) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	t_FluxEu flux;

	if (face.BCId.get() == t_FaceBCID::Fluid) {

		t_ConsVars csv_op = getCellCSV(iZone, face.pOppCell->Id);

		t_PrimVars pvl = csv_my.calcPrimVars();
		t_PrimVars pvr = csv_op.calcPrimVars();

		// TODO: verify
		// distance along normal between cell centers 
		double dx = abs((face.pOppCell->Center - face.pMyCell->Center).dot(face.Normal));
		double dt_dx = this->timeStep/dx;

		// rotate everything to local rf
		R.set(mat_rot_coefs);

		// debug
		//hsLogMessage("Face #%d:", iFace);
		//hsLogMessage(R.to_str().c_str());

		pvl.rotate(R);
		pvr.rotate(R);

		calcLaxWendroffFlux(pvl, pvr, dt_dx, flux);

		// rotate flux back
		R.set_inv(mat_rot_coefs);
		flux.rotate(R);

		// set flux for the face
		getFlux(iZone, iFace) = flux;

		return;


	}

	// TODO: now if not fluid, it's a BC

	t_BCKindEuler bc_kind = G_BCListEuler.getKind(face.BCId.get());

	if (bc_kind == t_BCKindEuler::Inflow) {

		t_ConsVars csv_face;
		csv_face.setValAtInf();

		R.set(mat_rot_coefs);
		csv_face.rotate(R);

		flux.calc(csv_face);

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;
		return;

	}

	if (bc_kind == t_BCKindEuler::Outflow) {

		t_ConsVars csv_face;

		csv_face = csv_my;

		R.set(mat_rot_coefs);

		csv_face.rotate(R);

		flux.calc(csv_face);

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;
		return;

	}

	if ((bc_kind == t_BCKindEuler::Wall) || (bc_kind == t_BCKindEuler::Sym)) {

		t_ConsVars csv_virt;

		R.set(mat_rot_coefs);

		csv_my.rotate(R);

		csv_virt = csv_my;
		// TODO: is this correct ? 
		csv_virt[1] *= -1.0;

		t_PrimVars pv_loc_my = csv_my.calcPrimVars();

		t_PrimVars pv_loc_virt = csv_virt.calcPrimVars();

		// distance between cell center and mirror cell center
		double dx = 2.0*abs((face.Center - face.pMyCell->Center).dot(face.Normal));

		double dt_dx = this->timeStep/dx;

		calcLaxWendroffFlux(pv_loc_my, pv_loc_virt, dt_dx, flux);

		// DEBUG, get normal velocity (normal to the wall)
		// if local rf this is u_ksi
		// IMPORTANT: instead of u_ksi, we now check (rho*u_ksi)
		// it is easy, TODO: compute u_ksi via flux vars
		double v_norm = flux[0];
		if (v_norm > G_State.ResidNormVeloWall) G_State.ResidNormVeloWall = v_norm;

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;

		return;

	}


	hsLogError(
		"t_DomainEuler::calcFaceFlux: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);



}


void calcLaxWendroffFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, const double dt_dx, t_FluxEu& flux) {

	t_ConsVars cvl, cvr;

	t_FluxEu fl, fr;

	calcCVFlux(pvl, cvl, fl);
	calcCVFlux(pvr, cvr, fr);

	// calculate LW average : 

	t_PrimVars pvm = 0.5 * (pvl + pvr);

	const double& u = pvm.getU();
	const double& v = pvm.getV();
	const double& w = pvm.getW();

	const double& Gamma = G_FlowModelParams.Gamma;

	const double K = Gamma - 1.0;

	const double q2 = u * u + v * v + w * w;
	// TODO: H: this depends on flow model
	const double H = pvm.calcH();

	t_SqMat<NConsVars> A;

	A[0][1] = 1;

	A[1][0] = 0.5 * K * q2 - u * u;
	A[1][1] = (3 - Gamma) * u;
	A[1][2] = (1 - Gamma) * v;
	A[1][3] = (1 - Gamma) * w;
	A[1][4] = Gamma - 1.0;

	A[2][0] = -u * v;
	A[2][1] = v;
	A[2][2] = u;

	A[3][0] = -u * w;
	A[3][1] = w;
	A[3][3] = u;

	A[4][0] = (0.5 * K * q2 - H) * u;
	A[4][1] = H + (1 - Gamma) * u * u;
	A[4][2] = (1 - Gamma) * u * v;
	A[4][3] = (1 - Gamma) * u * w;
	A[4][4] = Gamma * u;

	t_ConsVars dE_LW = dt_dx* (A * (fr - fl));
	// flux = 0.5 * (fl + fr - dt_dx * A * (fr - fl));
	for (int i = 0; i < NConsVars; i++)
		flux[i] = 0.5 * (fl[i] + fr[i] - dE_LW[i]);

};