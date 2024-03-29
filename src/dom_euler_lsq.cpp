#include "dom_euler_lsq.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "ghost_euler.h"

void t_DomEuLSQ::allocateFlowSolution() {

	t_DomEuBase::allocateFlowSolution();

	t_LSQData::allocateLSQData(this);

}

void t_DomEuLSQ::prepareBeforeTimeMarch() {

	t_LSQData::calcReconstData(G_GhostMngEu);

}

t_ConsVars t_DomEuLSQ::calcVirtCellCSV(int iZone, lint iFace) const{

	t_ConsVars csv_virt;

	const t_Zone& zne = Zones[iZone];

	const t_Face& face = zne.getFace(iFace);

	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	{
		TLogSyncGuard sg;
		if (face.isFluid()) 
			hsLogError("t_DomEuLSQ::calcVirtCellCSV: trying to calc virt cell csv for fluid face");
	}
		
	t_BCKindEuler bc_kind = G_BCListEuler.getKind(face.BCId.get());

	if (bc_kind == t_BCKindEuler::Inflow) {

		csv_virt.setValAtInf();
		return csv_virt;

	}

	if (bc_kind == t_BCKindEuler::Outflow) {

		csv_virt = csv_my;
		return csv_virt;

	}

	if ((bc_kind == t_BCKindEuler::Wall) || (bc_kind == t_BCKindEuler::Sym)) {

		t_MatRotN mat_rot_coefs;

		mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

		t_SqMat3 R;

		R.set(mat_rot_coefs);

		csv_my.rotate(R);

		csv_virt = csv_my;
		// TODO: is this correct ? 
		// flipping normal velocity
		csv_virt[1] *= -1.0;

		R.set_inv(mat_rot_coefs);

		csv_virt.rotate(R);
		return csv_virt;

	}


	hsLogError(
		"t_DomEuLSQ::calcVirtCellCSV: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);

	return csv_virt;

}

void t_DomEuLSQ::calcFaceFlux(int iZone, lint iFace) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	t_FluxEu flux;

	t_ConsVars csv_op;
	if (face.isFluid())
		csv_op = getCellCSV(iZone, face.pOppCell->Id);
	else
		csv_op = calcVirtCellCSV(iZone, iFace);

	const t_Cell& CellMy = *face.pMyCell;
	const t_Cell& CellOp = *face.pOppCell;

	t_Mat<NConsVars, 3> CellGradCSVMy;
	t_Vec<NConsVars> limMy;
	calcCellGradCSV(iZone, face.pMyCell->Id, CellGradCSVMy);
	// TODO: ghost cells must receive!
	// virt cells grads are zero
	t_Mat<NConsVars, 3> CellGradCSVOp;
	t_Vec<NConsVars> limOp({0,0,0,0,0});
	if (face.isFluid()) {
		calcCellGradCSV(iZone, face.pOppCell->Id, CellGradCSVOp);
		calcSlopeLimiters(iZone, face.pOppCell->Id, CellGradCSVOp);
	}

	// distances from cell centers to face center
	t_Vec3 drMy = face.Center - CellMy.Center;
	t_Vec3 drOp = face.Center - CellOp.Center;

	t_ConsVars dUMy = CellGradCSVMy * drMy;
	t_ConsVars dUOp = CellGradCSVOp * drOp;

	// apply limiters
	for (int k = 0; k < NConsVars; k++) {
		dUMy[k] *= limMy[k];
		dUOp[k] *= limOp[k];
	}

	t_ConsVars csv_l = csv_my + dUMy;
	t_ConsVars csv_r = csv_op + dUOp;

	t_PrimVars pvl = csv_l.calcPrimVars();
	t_PrimVars pvr = csv_r.calcPrimVars();

	// rotate everything to local rf
	R.set(mat_rot_coefs);

	pvl.rotate(R);
	pvr.rotate(R);

	calcRSFlux(pvl, pvr, flux);

	// rotate flux back
	R.set_inv(mat_rot_coefs);
	flux.rotate(R);

	// set flux for the face
	getFlux(iZone, iFace) = flux;

	return;

};