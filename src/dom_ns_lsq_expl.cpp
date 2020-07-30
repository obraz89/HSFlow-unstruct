#include "dom_ns_lsq_expl.h"

#include "flux_ns.h"

#include "rs_euler.h"

#include "bc_ns.h"

#include "ghost_ns.h"

void t_DomNSLSQ::allocateFlowSolution() {

	t_DomNSBase::allocateFlowSolution();

	t_LSQData::allocateLSQData(this);

}

void t_DomNSLSQ::prepareBeforeTimeMarch() {


	hsLogMessage("Calculating weights for vertices");

	calcCellWeightsForVertices();

	hsLogMessage("Calculating face grad reconstruction matrices");

	calcFaceGradMatrices();

	hsLogMessage("Calculating lsq reconstruction matrices");

	t_LSQData::calcReconstData(G_GhostMngNS);

}

// calculate csv in virtual cell depending on bc
// globar reference frame
t_ConsVars t_DomNSLSQ::calcVirtCellCSV(int iZone, lint iFace) const {

	t_ConsVars csv_virt;

	const t_Zone& zne = Zones[iZone];

	const t_Face& face = zne.getFace(iFace);

	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	{
		TLogSyncGuard sg;
		if (face.isFluid())
			hsLogError("t_DomEuLSQ::calcVirtCellCSV: trying to calc virt cell csv for fluid face");
	}

	int bc_id = face.BCId.get();

	t_BCKindNS bc_kind = G_BCListNS.getKind(bc_id);

	if (bc_kind == t_BCKindNS::InflowSup) {

		csv_virt.setValAtInf();
		return csv_virt;

	}

	if (bc_kind == t_BCKindNS::OutflowSup) {

		csv_virt = csv_my;
		return csv_virt;

	}

	if (bc_kind == t_BCKindNS::EulerWall ) {

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

	if (bc_kind == t_BCKindNS::WallNoSlip) {

		double Tw = G_BCListNS.getBC(bc_id)->get_settings_grp("").get_real_param("Tw");

		t_PrimVars pv_face = csv_my.calcPrimVars();

		// rho 
		pv_face.setR(calcRhoByPT(pv_face.getP(), Tw));
		// u,v,w
		pv_face.setUVW({0,0,0});
		// using dp/dn=0, P_face = P_my

		return pv_face.calcConsVars();

	}

	hsLogError(
		"t_DomNSLSQ::calcVirtCellCSV: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);

	return csv_virt;

}

void t_DomNSLSQ::calcFaceFluxEuler(int iZone, lint iFace, t_FluxEu& fluxEU) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	fluxEU.reset();

	t_ConsVars csv_op;
	if (face.isFluid())
		csv_op = getCellCSV(iZone, face.pOppCell->Id);
	else
		csv_op = calcVirtCellCSV(iZone, iFace);

	const t_Cell& CellMy = *face.pMyCell;
	const t_Cell& CellOp = *face.pOppCell;

	//compute inviscid flux
	
	t_Mat<NConsVars, 3> CellGradCSVMy;
	t_Vec<NConsVars> limMy;
	calcCellGradCSV(iZone, face.pMyCell->Id, CellGradCSVMy);
	// TODO: ghost cells must receive!
	// virt cells grads are zero
	t_Mat<NConsVars, 3> CellGradCSVOp;
	t_Vec<NConsVars> limOp({ 0,0,0,0,0 });
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

	calcRSFlux(pvl, pvr, fluxEU);

	// rotate flux back
	R.set_inv(mat_rot_coefs);
	fluxEU.rotate(R);

	return;

};

void t_DomNSLSQ::calcDataForFaceGradRUVWT(int iZone, lint iFace, t_VecConsVars& Umy,
	t_VecConsVars& Uop, t_Mat<NConsVars, MaxNumVertsInFace>& UVerts) const{
	
	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	int bcid = face.BCId.get();

	t_PrimVars pv_my = getCellCSV(iZone, face.pMyCell->Id).calcPrimVars();

	t_PrimVars pv_op;

	if (face.isFluid())
		pv_op = getCellCSV(iZone, face.pOppCell->Id).calcPrimVars();
	else
		pv_op = calcVirtCellCSV(iZone, iFace).calcPrimVars();

	Umy = pv_my.calcRUVWT();

	Uop = pv_op.calcRUVWT();

	t_PrimVars pv_vert;

	for (int ivert = 0; ivert < face.NVerts; ivert++) {

		pv_vert = getVertCSV(iZone, face.pVerts[ivert]->Id).calcPrimVars();
		UVerts.setCol(ivert, pv_vert.calcRUVWT());

	}

	return;

};

void t_DomNSLSQ::calcFaceFluxVisc(int iZone, lint iFace, t_VecConsVars& fluxVisc) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);


	t_VecConsVars Umy, Uop;
	t_Mat<NConsVars, MaxNumVertsInFace> UVerts;
	calcDataForFaceGradRUVWT(iZone, iFace, Umy, Uop, UVerts);

	t_Mat<3, NConsVars> FaceGradRUVWT;
	face.ComputeFaceGrad<5>(Umy, Uop, UVerts, FaceGradRUVWT);


	// IMPORTANT TODO: prim values at face must be computed as averaged 
	// reconstructed values at face center: 
	// pv_face = 0.5*(pv_left_lsq + pv_right_lsq)
	
	// for now just using cell center value!
	t_PrimVars pv_face = getCellCSV(iZone, face.pMyCell->Id).calcPrimVars();

	calcNSViscFlux(face.Normal, pv_face, FaceGradRUVWT, fluxVisc);

}

void t_DomNSLSQ::calcFaceFlux(int iZone, lint iFace) {

	t_VecConsVars fluxTot;

	t_FluxEu fluxEu;

	calcFaceFluxEuler(iZone, iFace, fluxEu);

	t_VecConsVars fluxVisc;

	calcFaceFluxVisc(iZone, iFace, fluxVisc);

	// set flux for the face

	fluxTot = fluxEu; //+ fluxVisc;

	getFlux(iZone, iFace) = fluxTot;

}