#include "dom_euler_lsq.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

void t_DomEuLSQ::allocateFlowSolution() {

	t_DomEuBase::allocateFlowSolution();

	ZonesRecData = new t_ZoneReconstData[nZones];

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneReconstData& fdata = ZonesRecData[i];

		lint nCellsTot = zne.getnCellsTot();

		fdata.ReconstData = new t_ReconstDataLSQ[nCellsTot];


	}
}

t_DomEuLSQ::~t_DomEuLSQ() {

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneReconstData& fdata = ZonesRecData[i];

		delete[] fdata.ReconstData;

	}

	delete[] ZonesRecData;

}

static void setDmMat(const t_Vec3& dr, t_SqMat3& dM) {

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dM[i][j] = dr[i] * dr[j];

}

void t_DomEuLSQ::calcReconstDataLSQ(int iZone, lint iCell) {

	t_SqMat3 M, dM;

	t_Vec3 dr;

	const t_Zone& Zne = Zones[iZone];

	const t_Cell& Cell = Zne.getCell(iCell);

	t_ReconstDataLSQ& ReconstData = getReconstData(iZone, iCell);

	const t_Vec3& xc = Cell.Center;

	// first find max distance btw cells, this is a scale factor r

	const int NNeig = Cell.NCellsNeig();

	double r = 0.0;

	for (int j = 0; j < NNeig; j++) {
		dr = Cell.pCellsNeig[j]->Center - xc;
		double r_cur = dr.norm();
		if (r_cur > r) r = r_cur;
	}

	ReconstData.r = r;

	double r_inv = 1.0 / r;

	for (int j = 0; j < Cell.NCellsNeig(); j++) {

		dr = Cell.pCellsNeig[j]->Center - xc;

		dr *= r_inv;

		setDmMat(dr, dM);

		M.add(dM);

	}

	ReconstData.MInvRR = M.CalcInv();

};

void t_DomEuLSQ::calcReconstData() {

	for (int iZone = 0; iZone < nZones; iZone++) {

		const t_Zone& zne = Zones[iZone];

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			calcReconstDataLSQ(iZone, iCell);

		}

	}

}

void t_DomEuLSQ::calcFaceFlux(int iZone, lint iFace) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	t_Flux flux;

	if (face.BCId.get() == t_FaceBCID::Fluid) {

		t_ConsVars csv_op = getCellCSV(iZone, face.pOppCell->Id);

		t_PrimVars pvl = csv_my.calcPrimVars();
		t_PrimVars pvr = csv_op.calcPrimVars();

		// rotate everything to local rf
		R.set(mat_rot_coefs);

		// debug
		//hsLogMessage("Face #%d:", iFace);
		//hsLogMessage(R.to_str().c_str());

		pvl.rotate(R);
		pvr.rotate(R);

		calcRusanovFlux(pvl, pvr, flux);

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
		// TODOL is this correct ? 
		csv_virt[1] *= -1.0;

		t_PrimVars pv_loc_my = csv_my.calcPrimVars();

		t_PrimVars pv_loc_virt = csv_virt.calcPrimVars();

		calcRusanovFlux(pv_loc_my, pv_loc_virt, flux);

		// DEBUG, get normal velocity (normal to the wall)
		// if local rf this is u_ksi
		// IMPORTANT: instead of u_ksi, we now check (rho*u_ksi)
		// it is easy, TODO: compute u_ksi via flux vars
		double v_norm = flux[0];
		if (v_norm > G_State.ResidNormVeloWall) G_State.ResidNormVeloWall = v_norm;

		// TODO: is this correct
		//flux[1] = 0.0;

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;

		return;

	}


	hsLogError(
		"t_DomainEuler::calcFaceFlux: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);

};