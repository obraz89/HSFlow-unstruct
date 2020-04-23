#include "reconstruct_lsq.h"

void setDmMat(const t_Vec3& dr, t_SqMat3& dM) {

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dM[i][j] = dr[i] * dr[j];

}

t_ReconstDataLSQ calcReconstDataLSQ(const t_Cell& Cell) {

	t_ReconstDataLSQ ReconstrData;

	t_SqMat3 M, dM;

	t_Vec3 dr;

	const t_Vec3& xc = Cell.Center;

	// first find max distance btw cells, this is a scale factor r

	const int NNeig = Cell.NCellsNeig();

	double r = 0.0;

	for (int j = 0; j < NNeig; j++) {
		dr = Cell.pCellsNeig[j]->Center - xc;
		double r_cur = dr.norm();
		if (r_cur > r) r = r_cur;
	}

	ReconstrData.r = r;

	double r_inv = 1.0 / r;
	
	for (int j = 0; j < Cell.NCellsNeig(); j++) {

		dr = Cell.pCellsNeig[j]->Center - xc;

		dr *= r_inv;

		setDmMat(dr, dM);

		M.add(dM);

	}

	ReconstrData.MInvRR = M.CalcInv();

	return ReconstrData;

};