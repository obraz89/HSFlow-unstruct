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

	// first find avg distance btw cells, this is a scale factor r
	// TODO: best scaling? avg, min, max?

	const int NNeig = Cell.NCellsNeig();

	double r = 0.0;

	for (int j = 0; j < Cell.NFaces; j++) {

		if (Cell.pCellsNeig[j] != nullptr) {

			dr = Cell.pCellsNeig[j]->Center - xc;
			r += dr.norm();

		}
	}

	r /= double(NNeig);

	ReconstData.r = r;

	double r_inv = 1.0 / r;

	for (int j = 0; j < Cell.NFaces; j++) {

		if (Cell.pCellsNeig[j] != nullptr) {

			dr = Cell.pCellsNeig[j]->Center - xc;

			dr *= r_inv;

			setDmMat(dr, dM);

			M.add(dM);

		}

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
// we have precomputed scaled inverse Mc matrix
// grad(U_i) = McInv*Summ((rd-rc)/r*(Uc_i - Ud_i)/r)
// i is component of prim or consv vars to reconstruct i=0...NConsVars-1
// r is scale introduced to make all multipliers in formula O(1)
void t_DomEuLSQ::calcCellGradPrimVars(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGradPV) {

	const t_ConsVars& csv_c = getCellCSV(iZone, iCell);

	const t_Cell& Cell = Zones[iZone].getCell(iCell);

	const t_ReconstDataLSQ& RecData = getReconstData(iZone, iCell);

	t_PrimVars pvs_c = csv_c.calcPrimVars();

	CellGradPV.reset();
	
	t_Vec3 grad_cur, dr;

	double du;

	double r_inv = 1.0 / RecData.r;

	t_PrimVars pvs_n;

	for (int j = 0; j < Cell.NFaces; j++) {

		if (Cell.pCellsNeig[j] != nullptr) {

			const t_Cell& CellNeig = *Cell.pCellsNeig[j];

			pvs_n = getCellCSV(iZone, CellNeig.Id).calcPrimVars();

			dr = r_inv * (CellNeig.Center - Cell.Center);

			for (int i = 0; i < NConsVars; i++) {

				du = r_inv * (pvs_c[i] - pvs_n[i]);

				grad_cur = RecData.MInvRR * dr;
				grad_cur *= du;

				for (int k = 0; k < 3; k++)
					CellGradPV[i][k] += grad_cur[k];
			}

		}

	}

};

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

		const t_Cell& CellMy = Zones[iZone].getCell(face.pMyCell->Id);
		const t_Cell& CellOp = Zones[iZone].getCell(face.pOppCell->Id);

		// prim vars at cell centers
		t_PrimVars pvl_c = csv_my.calcPrimVars();
		t_PrimVars pvr_c = csv_op.calcPrimVars();

		t_Mat<NConsVars, 3> CellGradPVMy;
		calcCellGradPrimVars(iZone, face.pMyCell->Id, CellGradPVMy);
		// TODO:: compute only for real cells, ghost cells must receive!
		t_Mat<NConsVars, 3> CellGradPVOp;
		calcCellGradPrimVars(iZone, face.pOppCell->Id, CellGradPVOp);

		// distances from cell centers to face center
		t_Vec3 drMy = face.Center - CellMy.Center;
		t_Vec3 drOp = face.Center - CellOp.Center;

		t_PrimVars dUMy = CellGradPVMy * drMy;
		t_PrimVars dUOp = CellGradPVOp * drOp;

		t_PrimVars pvl = pvl_c + dUMy;
		t_PrimVars pvr = pvr_c + dUOp;

		// rotate everything to local rf
		R.set(mat_rot_coefs);

		// debug
		//hsLogMessage("Face #%d:", iFace);
		//hsLogMessage(R.to_str().c_str());

		pvl.rotate(R);
		pvr.rotate(R);

		calcRSFlux(pvl, pvr, flux);

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

		calcRSFlux(pv_loc_my, pv_loc_virt, flux);

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

};