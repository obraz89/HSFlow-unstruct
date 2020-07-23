#include "lsq_uvwpt.h"

void t_LSQData::allocateLSQData(t_Dom5* a_dom) {

	Dom = a_dom;

	ZonesRecData = new t_ZoneReconstData[Dom->nZones];
	ZonesVirtCells = new t_ZoneVirtCells[Dom->nZones];

	for (int i = Dom->iZneMPIs; i <= Dom->iZneMPIe; i++) {

		t_Zone& zne = Dom->Zones[i];

		ZonesRecData[i].ReconstData = new t_ReconstDataLSQ[zne.getnCellsTot()];

		ZonesVirtCells[i].VirtCells = new t_Cell[zne.getNFacesBC()];

	}

	initializeVirtCells();

}
// attach virtual cells to real cells at BC faces
// after this all real cells should have full set of neighbors;
// virtual cells has only center and cell center values
void t_LSQData::initializeVirtCells() {
	for (int i = Dom->iZneMPIs; i <= Dom->iZneMPIe; i++) {

		t_Zone& zne = Dom->Zones[i];

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

t_LSQData::~t_LSQData() {

	for (int i = Dom->iZneMPIs; i <= Dom->iZneMPIe; i++) {

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

void t_LSQData::calcReconstDataLSQ(int iZone, lint iCell) {

	t_SqMat3 M, dM;

	t_Vec3 dr;

	const t_Zone& Zne = Dom->Zones[iZone];

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

void t_LSQData::calcReconstData(t_GhostMng5& ghost_mng) {

	for (int iZone = Dom->iZneMPIs; iZone <= Dom->iZneMPIe; iZone++) {

		const t_Zone& zne = Dom->Zones[iZone];

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			calcReconstDataLSQ(iZone, iCell);

		}

	}

	ghost_mng.exchangeReconstData();

}
// we have precomputed scaled inverse Mc matrix;
// grad(U_i) = McInv*Summ((rd-rc)/r*(Ud_i - Uc_i)/r)
// i is component of consv vars to reconstruct i=0...NConsVars-1
// r is scale introduced to make all multipliers in formula O(1)
void t_LSQData::calcCellGradCSV(int iZone, lint iCell, t_Mat<NConsVars, 3>& CellGradCSV) {

	const t_ConsVars& csv_c = Dom->getCellCSV(iZone, iCell);

	const t_Cell& Cell = Dom->Zones[iZone].getCell(iCell);

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

static double calcLimiterMinmod(double r) {
	return fmin(fabs(r), 1.0);
}

t_Vec<NConsVars> t_LSQData::calcSlopeLimiters(
	int iZone, lint iCell, const t_Mat<NConsVars, 3>& CellGradCSV) const {

	t_Vec<NConsVars> limiters({ 1,1,1,1,1 });

	const t_Cell& Cell = Dom->Zones[iZone].getCell(iCell);

	const t_ConsVars& csvMy = Dom->getCellCSV(iZone, iCell);

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

// get csv of neighbor cell
// if it is real cell, just read the values
// if it is virtual cell, compute csv
t_ConsVars t_LSQData::getNeigCellCSV(int iZone, lint iCell, int indFace) const {

	const t_Cell& Cell = Dom->Zones[iZone].getCell(iCell);

	const t_Face& Face = Cell.getFace(indFace);

	t_ConsVars csv;

	if (Face.isFluid())
		csv = Dom->getCellCSV(iZone, Cell.pCellsNeig[indFace]->Id);
	else
		csv = calcVirtCellCSV(iZone, Face.Id);

	return csv;
};