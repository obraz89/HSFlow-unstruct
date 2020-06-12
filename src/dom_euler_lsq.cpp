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
// i is component of consv vars to reconstruct i=0...NConsVars-1
// r is scale introduced to make all multipliers in formula O(1)
void t_DomEuLSQ::calcCellGradCSV(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGradCSV) {

	const t_ConsVars& csv_c = getCellCSV(iZone, iCell);

	const t_Cell& Cell = Zones[iZone].getCell(iCell);

	const t_ReconstDataLSQ& RecData = getReconstData(iZone, iCell);

	CellGradCSV.reset();
	
	t_Vec3 grad_cur, dr;

	double du;

	double r_inv = 1.0 / RecData.r;

	t_ConsVars csv_n;

	// iterate over neighbors, including virtual cells
	for (int j = 0; j < Cell.NFaces; j++) {

		const t_Cell& CellNeig = *Cell.pCellsNeig[j];

		dr = r_inv * (CellNeig.Center - Cell.Center);

		csv_n = getNeigCellCSV(iZone, iCell, j);

		for (int i = 0; i < NConsVars; i++) {

			du = r_inv * (csv_n[i] - csv_c[i]);

			grad_cur = RecData.MInvRR * dr;
			grad_cur *= du;

			for (int k = 0; k < 3; k++)
				CellGradCSV[i][k] += grad_cur[k];
		}

	}

};

double calcLimiterMinmod(double r) {
	return fmin(fabs(r), 1.0);
}

t_Vec<NConsVars> t_DomEuLSQ::calcSlopeLimiters(
	int iZone, lint iCell, const t_Mat<NConsVars, 3>& CellGradCSV) const {

	t_Vec<NConsVars> limiters({1,1,1,1,1});

	const t_Cell& Cell = Zones[iZone].getCell(iCell);

	const t_ConsVars& csvMy = getCellCSV(iZone, iCell);

	t_ConsVars csvMin = csvMy;
	t_ConsVars csvMax = csvMy;

	t_ConsVars csvNeig;

	// iterate over neighbors to find min max values
	for (int j = 0; j < Cell.NFaces; j++) {

		csvNeig = getNeigCellCSV(iZone, iCell, j);

		for (int k = 0; k < NConsVars; k++) {
			if (csvNeig[k] > csvMax[k]) csvMax[k] = csvNeig[k];
			if (csvNeig[k] < csvMin[k]) csvMin[k] = csvNeig[k];
		}

	}

	t_Vec3 dr;

	t_ConsVars csvVert;

	{
		double r;
		double lim; 

		// if difference in fluid vars is less then tol 
		// limiter is not applied
		// TODO: this is an empirical constant
		const double LIM_TOL = 1.0e-03;

		// iterate over vertices, 
		// compute limiter for each vertex
		double Uc;
		double Up;
		double Umax, Umin;
		for (int nv = 0; nv < Cell.NVerts; nv++) {

			dr = Cell.getVert(nv).xyz - Cell.Center;

			csvVert = csvMy + CellGradCSV * dr;

			for (int k = 0; k < NConsVars; k++) {

				Uc = csvMy[k];
				Up = csvVert[k];

				Umax = csvMax[k];
				Umin = csvMin[k];

				if (fabs(Up - Uc) > LIM_TOL) {
					if (Up > Uc)
						r = (Umax - Uc) / (Up - Uc);
					else
						r = (Umin - Uc) / (Up - Uc);

					lim = calcLimiterMinmod(r);
				}
				else {
					lim = 1.0;
				}

				if (lim < limiters[k]) limiters[k] = lim;

			}

		}

	}

	return limiters;

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

// get csv of neighbor cell
// if it is real cell, just read the values
// if it is virtual cell, compute csv
t_ConsVars t_DomEuLSQ::getNeigCellCSV(int iZone, lint iCell, int indFace) const {

	const t_Cell& Cell = Zones[iZone].getCell(iCell);

	const t_Face& Face = Cell.getFace(indFace);

	t_ConsVars csv;

	if (Face.isFluid())
		csv = getCellCSV(iZone, Cell.pCellsNeig[indFace]->Id);
	else
		csv = calcVirtCellCSV(iZone, Face.Id);

	return csv;
};

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