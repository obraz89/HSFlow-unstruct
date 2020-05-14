#include "dom_euler_lsq.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "ghost_euler.h"

void t_DomEuLSQ::allocateFlowSolution() {

	t_DomEuBase::allocateFlowSolution();

	ZonesRecData = new t_ZoneReconstData[nZones];
	ZonesVirtCells = new t_ZoneVirtCells[nZones];

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];

		ZonesRecData[i].ReconstData = new t_ReconstDataLSQ[zne.getnCellsTot()];
		
		ZonesVirtCells[i].VirtCells = new t_Cell[zne.getNFacesBC()];

	}

	initializeVirtCells();

}
// attach virtual cells to real cells at BC faces
// after this all real cells should have full set of neighbors;
// virtual cells has only center and cell center values
void t_DomEuLSQ::initializeVirtCells() {
	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];

		t_Cell* VirtCells = ZonesVirtCells[i].VirtCells;

		int iCellVirt = 0;

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			t_Cell& cell = zne.getCell(iCell);

			for (int indFace = 0; indFace < cell.NFaces; indFace++) {
				const t_Face& face = cell.getFace(indFace);
				if (face.isFluid()) continue;
				// to update face we need access via zone
				t_Face& face2change = zne.getFace(face.Id);

				t_Cell& VirtCell = VirtCells[iCellVirt++];

				cell.pCellsNeig[indFace] = &VirtCell;
				face2change.pOppCell = &VirtCell;

				// compute virtual cell center
				// flipping real cell center over face

				const t_Vec3& n = face.Normal;

				double dst = 2.0 * n.dot(face.Center - cell.Center);

				VirtCell.Center = cell.Center + dst * n;

			}
		}

	}
}

t_DomEuLSQ::~t_DomEuLSQ() {

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		delete[] ZonesRecData[i].ReconstData;
		delete[] ZonesVirtCells[i].VirtCells;

	}

	delete[] ZonesRecData;
	delete[] ZonesVirtCells;

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

	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

		const t_Zone& zne = Zones[iZone];

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			calcReconstDataLSQ(iZone, iCell);

		}

	}

	G_GhostMngEu.exchangeReconstData();

}
// we have precomputed scaled inverse Mc matrix;
// grad(U_i) = McInv*Summ((rd-rc)/r*(Ud_i - Uc_i)/r)
// i is component of prim or consv vars to reconstruct i=0...NConsVars-1
// r is scale introduced to make all multipliers in formula O(1)
void t_DomEuLSQ::calcCellGradPrimVars(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGradPV) {

	const t_ConsVars& csv_c = getCellCSV(iZone, iCell);

	const t_Cell& Cell = Zones[iZone].getCell(iCell);

	const t_ReconstDataLSQ& RecData = getReconstData(iZone, iCell);

	t_PrimVars pvs_c = csv_c.calcPrimVars();

	CellGradPV.reset();

	return;
	
	t_Vec3 grad_cur, dr;

	double du;

	double r_inv = 1.0 / RecData.r;

	t_PrimVars pvs_n;

	// iterate over neighbors, including virtual cells
	for (int j = 0; j < Cell.NFaces; j++) {

		const t_Face& face = Cell.getFace(j);

		const t_Cell& CellNeig = *Cell.pCellsNeig[j];

		if (face.isFluid())
			pvs_n = getCellCSV(iZone, CellNeig.Id).calcPrimVars();
		else
			pvs_n = calcVirtCellCSV(iZone, face.Id).calcPrimVars();

		dr = r_inv * (CellNeig.Center - Cell.Center);

		for (int i = 0; i < NConsVars; i++) {

			du = r_inv * (pvs_n[i] - pvs_c[i]);

			grad_cur = RecData.MInvRR * dr;
			grad_cur *= du;

			for (int k = 0; k < 3; k++)
				CellGradPV[i][k] += grad_cur[k];
		}

	}

};

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

	t_Flux flux;

	t_ConsVars csv_op;
	if (face.isFluid())
		csv_op = getCellCSV(iZone, face.pOppCell->Id);
	else
		csv_op = calcVirtCellCSV(iZone, iFace);

	const t_Cell& CellMy = *face.pMyCell;
	const t_Cell& CellOp = *face.pOppCell;

	// prim vars at cell centers
	t_PrimVars pvl_c = csv_my.calcPrimVars();
	t_PrimVars pvr_c = csv_op.calcPrimVars();

	t_Mat<NConsVars, 3> CellGradPVMy;
	calcCellGradPrimVars(iZone, face.pMyCell->Id, CellGradPVMy);
	// TODO: ghost cells must receive!
	// virt cells grads are zero
	t_Mat<NConsVars, 3> CellGradPVOp;
	if (face.isFluid())
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