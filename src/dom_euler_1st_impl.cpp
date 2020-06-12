#include "dom_euler_1st_impl.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "ghost_euler.h"

int t_DomEu1stImpl::getLocInd(const int iZone, const int iCell) const{

#ifdef _DEBUG
	if (!Zones[iZone].isRealCell(iCell))
		hsLogError("t_DomEu1stImpl::getLocInd: trying to calc loc index for ghost cell (cell must be real!)");
#endif

	int idZoneOffset = 0;

	for (int i = iZneMPIs; i < iZone; i++)
		idZoneOffset += Zones[i].getnCellsReal();

	return idZoneOffset + iCell;


}

int t_DomEu1stImpl::getGlobInd(const int iZone, const int iCell) const{

#ifdef _DEBUG
	if (iZone<iZneMPIs || iZone>iZneMPIe)
		hsLogError("t_DomEu1stImpl::getGlobInd: trying to calculate global index for cell that is not in my zones");
#endif

	const t_Zone& zne = Zones[iZone];

	if (zne.isRealCell(iCell)) {

		int idZoneOffset = 0;

		for (int i = 0; i < iZone; i++)
			idZoneOffset += Zones[i].getnCellsReal();

		return idZoneOffset + iCell;

	}
	else {
		int iZoneDnr;
		int iCellDnr;
		G_GhostMngEu.getDonor(iZone, iCell, iZoneDnr, iCellDnr);

		int idZoneOffset = 0;

		for (int i = 0; i < iZoneDnr; i++)
			idZoneOffset += Zones[i].getnCellsReal();

		return idZoneOffset + iCellDnr;

	}

};

bool t_DomEu1stImpl::isMyCell(int iZone, int iCell) const{

	if (iZone<iZneMPIs || iZone>iZneMPIe)
		hsLogError("t_DomEu1stImpl::isMyCell: trying to check cell that is not in my zones");

	const t_Zone& Zne = Zones[iZone];

	// check if cell is real
	if (Zne.isRealCell(iCell))
		return true;
	// check if ghost but donor also belongs to me
	int iZoneDnr;
	int iCellDnr;
	G_GhostMngEu.getDonor(iZone, iCell, iZoneDnr, iCellDnr);

	if (iZneMPIs <= iZoneDnr && iZoneDnr <= iZneMPIe)
		return true;

	// cell is ghost and donor is not mine
	return false;

}

void t_DomEu1stImpl::allocateFlowSolution() {

	t_DomEuBase::allocateFlowSolution();

	ZonesLambdaCD = new t_ZoneLambdasCD[nZonesMy()];

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {
		const t_Zone& zne = Zones[i];
		int nfaces = zne.getNFaces();

		ZonesLambdaCD[i].Lambdas = new double[nfaces];

		double* lambdas = ZonesLambdaCD[i].Lambdas;

		for (int j = 0; j < nfaces; j++)
			lambdas[j] = 0.0;
	}

	// get dimensions of global unknown vector du
	// and starting offset for each rank
	int NCellsMy = 0;

	for (int i = iZneMPIs; i <= iZneMPIe; i++)
			NCellsMy += Zones[i].getnCellsReal();

	ctxKSP.dimMy = NCellsMy*NConsVars;

	// dimGlob is a sum of dimMy for all ranks 
	MPI_Allreduce(&ctxKSP.dimMy, &ctxKSP.dimGlob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// number of non-zero blocks in diag block-row (d_nnz) and non-diagonal block-row (o_nnz) 
	// diag submatrix is for all my non-ghost cells
	// non-diagonal blocks occur when cell neighbor is ghost
	int* d_nnz = new int[NCellsMy];
	int* o_nnz = new int[NCellsMy];

	//debug
	int maxNDiagBlocks = 0;
	int maxNNonDiagBlocks = 0;

	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {
		const t_Zone& Zne = Zones[iZone];
		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {
			int iBlockRowLocal = getLocInd(iZone, iCell);
			const t_Cell& cell = Zne.getCell(iCell);
			int NDiagBlocks = 1;
			int NNonDiagBlocks = 0;
			for (int k = 0; k < MaxNumFacesInCell; k++) {
				if (cell.pCellsNeig[k] != nullptr) {
					const t_Cell& cell_op = *cell.pCellsNeig[k];
					if (isMyCell(iZone, cell_op.Id))
						NDiagBlocks++;
					else
						NNonDiagBlocks++;
				}
			}
			// number of non-diag blocks
			d_nnz[iBlockRowLocal] = NDiagBlocks;
			o_nnz[iBlockRowLocal] = NNonDiagBlocks;

			//debug
			if (NDiagBlocks > maxNDiagBlocks) maxNDiagBlocks = NDiagBlocks;
			if (NNonDiagBlocks > maxNNonDiagBlocks) maxNNonDiagBlocks = NNonDiagBlocks;
		}
	}

	hsLogMessage("MaxDiagBlocks=%d, MaxNonDiagBlocks=%d", maxNDiagBlocks, maxNNonDiagBlocks);

	PetscErrorCode perr;
	// create vectors
	// create rhs
	perr = VecCreateMPI(PETSC_COMM_WORLD, ctxKSP.dimMy, PETSC_DETERMINE, &ctxKSP.b);
	if (perr) hsLogMessage("Error:Failed to init PETSc rhs vector");

	// x - store exactly as rhs as we need to copy data from solution later
	perr = VecDuplicate(ctxKSP.b, &ctxKSP.x);
	if (perr) hsLogMessage("Error:Failed to init PETSc solution vector");

	// allocate matrix
	perr = MatCreateBAIJ(PETSC_COMM_WORLD, NConsVars,
		ctxKSP.dimMy, ctxKSP.dimGlob,
		PETSC_DETERMINE, ctxKSP.dimGlob,
		0, d_nnz, 0, o_nnz, &ctxKSP.A);
	if (perr) hsLogMessage("Error:Failed to init PETSc matrix");

	MatSetOption(ctxKSP.A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

	perr = KSPCreate(PETSC_COMM_WORLD, &ctxKSP.ksp);
	// deflated GMRES, the adaptive strategy allows to switch to the deflated GMRES when the stagnation occurs
	//perr = KSPSetType(ctxKSP.ksp, KSPDGMRES); 

	KSPSetOperators(ctxKSP.ksp, ctxKSP.A, ctxKSP.A);
	// ???
	//KSPSetTolerances(ctxKSP.ksp, 1.e-5, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetTolerances(ctxKSP.ksp, PETSC_DEFAULT, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT);

	//KSPSetFromOptions(ctxKSP.ksp);
	
	delete[] d_nnz, o_nnz;
}
// explicit face flux
void t_DomEu1stImpl::calcFaceFlux(int iZone, lint iFace) {

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

		// rotate everything to local rf
		R.set(mat_rot_coefs);

		// debug
		//hsLogMessage("Face #%d:", iFace);
		//hsLogMessage(R.to_str().c_str());

		pvl.rotate(R);
		pvr.rotate(R);

		calcRSFlux(pvl, pvr, flux);
		ZonesLambdaCD[iZone].Lambdas[face.Id] = 
			calcWaveSpeedDavisEstimMaxAbs(pvl,pvr);

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

		t_PrimVars pv = csv_face.calcPrimVars();
		// assuming that virtual cell pv is equal to cell pv
		ZonesLambdaCD[iZone].Lambdas[face.Id] =
			calcWaveSpeedDavisEstimMaxAbs(pv, pv);

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

		ZonesLambdaCD[iZone].Lambdas[face.Id] =
			calcWaveSpeedDavisEstimMaxAbs(pv_loc_my, pv_loc_virt);

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;

		return;

	}

	hsLogError(
		"t_DomainEuler::calcFaceFlux: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);

}

inline int plainInd4Blk(int i, int j) {
	return i * NConsVars + j;
}

void t_DomEu1stImpl::makeTimeStep() {

	double dt_local = calcDt();

	double dt;

	MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	double dU_max = 0.0;

	hsLogMessage("Computed dt=%lf", dt);

	G_State.ResidTot = 0.0;

	G_State.ResidNormVeloWall = 0.0;

	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

		t_Zone& zne = Zones[iZone];

		for (int iFace = 0; iFace < zne.getNFaces(); iFace++) {
			calcFaceFlux(iZone, iFace);
		}

		// constructing rhs (du from explicit scheme)
		{
			t_ConsVars dU_rhs;

			double buf[NConsVars];
			int idx[NConsVars];

			for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

				int idGlob = getGlobInd(iZone, iCell);
				int iRowBase = NConsVars * idGlob;

				dU_rhs.reset();

				t_Cell& cell = zne.getCell(iCell);

				for (int j = 0; j < cell.NFaces; j++) {

					const t_Face& face = cell.getFace(j);

					// U[n+1] = U[n] - summ(flux_j), so we summ fluxes multiplied by minus 1
					double coef = cell.isMyFace(j) ? -1.0 : 1.0;

					coef *= dt * face.Area / cell.Volume;

					const t_FluxEu& flux = getFlux(iZone, face.Id);

					dU_rhs += coef * flux;

				}

				// insert values in my portion of rhs vector
				// TODO:VecSetValuesBlocked
				for (int k = 0; k < NConsVars; k++) {
					buf[k] = dU_rhs[k];
					idx[k] = iRowBase + k;
				}

				VecSetValues(ctxKSP.b, NConsVars, idx, buf, INSERT_VALUES);

			}

		
		}

		// inserting values into The Matrix 
		{
			// reset non-zeros

			// one block per insertion
			int idxm;
			int idxn;
			double buf[NConsVars*NConsVars];

			t_MatRotN mat_rot_coefs;

			t_SqMat3 R, R_inv;

			t_SqMat<NConsVars> Jac_c, Jac_d;

			t_SqMat<NConsVars> Jac_c_glob, Jac_d_glob;

			t_SqMat<NConsVars> R_inv_infl, R_infl;

			for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

				const t_Cell& cell = zne.getCell(iCell);

				const int idGlob = getGlobInd(iZone, iCell);
				// diagonal unity block
				{
					for (int i = 0; i < NConsVars; i++) 
						for (int j = 0; j < NConsVars; j++)
							buf[plainInd4Blk(i, j)] = (i==j) ? 1.0 : 0.0;
				}
				
				idxm = idGlob;
				idxn = idGlob;
				int* pidxm = &idxm;
				int val = pidxm[0];
				MatSetValuesBlocked(ctxKSP.A, 1, &idxm, 1, &idxn, buf, ADD_VALUES);

				for (int j = 0; j < cell.NFaces; j++) {

					const t_Face& face = cell.getFace(j);

					t_Vec3 Normal = cell.getFaceNormalOutward(j);

					mat_rot_coefs.calc_rot_angles_by_N(Normal);

					R.set(mat_rot_coefs);
					R_inv.set_inv(mat_rot_coefs);

					t_ConsVars::inflateRotMat(mat_rot_coefs, R_infl);
					t_ConsVars::inflateRotMatInv(mat_rot_coefs, R_inv_infl);

					double coef = dt * (face.Area / cell.Volume);

					if (face.isFluid()) {

						t_PrimVars pvc, pvd;

						pvc = getCellCSV(iZone, iCell).calcPrimVars();

						const t_Cell& cell_opp = *cell.pCellsNeig[j];

						pvd = getCellCSV(iZone, cell_opp.Id).calcPrimVars();

						double LamAbs = ZonesLambdaCD[iZone].Lambdas[face.Id];

						// rotate to local rf

						pvc.rotate(R);
						pvd.rotate(R);

						pvc.calcJac(Jac_c);
						pvd.calcJac(Jac_d);

						// rotate matrices back
						// Jac_glob = R_inv*Jac_loc*R
						Jac_c_glob = R_inv_infl * Jac_c * R_infl;
						Jac_d_glob = R_inv_infl * Jac_d * R_infl;

						// diagonal block for c
						{
							double LamAdd;
							for (int i = 0; i < NConsVars; i++) {
								for (int j = 0; j < NConsVars; j++) {
									LamAdd = (i == j) ? LamAbs : 0.0;
									buf[plainInd4Blk(i, j)] = 
										0.5 * coef * (Jac_c_glob[i][j] + LamAdd);
								}
							}
							idxm = idGlob;
							idxn = idGlob;
							MatSetValuesBlocked(ctxKSP.A, 1, &idxm, 1, &idxn, buf, ADD_VALUES);
						}
						// offset block for d 
						{
							double LamAdd;
							for (int i = 0; i < NConsVars; i++) {
								for (int j = 0; j < NConsVars; j++) {
									LamAdd = (i == j) ? -1.0 * LamAbs : 0.0;
									buf[plainInd4Blk(i, j)] =
										0.5 * coef * (Jac_d_glob[i][j] + LamAdd);
								}
							}
							idxm = idGlob;
							idxn = getGlobInd(iZone, cell_opp.Id);
							MatSetValuesBlocked(ctxKSP.A, 1, &idxm, 1, &idxn, buf, ADD_VALUES);
						}

						continue;
					}

					t_BCKindEuler bc_kind = G_BCListEuler.getKind(face.BCId.get());

					if (bc_kind == t_BCKindEuler::Inflow) {
						// nothing to add, Fc=const, no jac additions
						continue;
					}

					if (bc_kind == t_BCKindEuler::Outflow ||
						bc_kind == t_BCKindEuler::Sym ||
						bc_kind == t_BCKindEuler::Wall) {
						// calculate diagonal addition, 
						// opposite cell is virtual, no non-diag addition
						t_PrimVars pvc;

						pvc = getCellCSV(iZone, iCell).calcPrimVars();

						double LamAbs = ZonesLambdaCD[iZone].Lambdas[face.Id];

						pvc.rotate(R);

						pvc.calcJac(Jac_c);

						// rotate Jac back
						Jac_c_glob = R_inv_infl * Jac_c * R_infl;

						// diagonal block for c
						{
							double LamAdd;
							for (int i = 0; i < NConsVars; i++) {
								for (int j = 0; j < NConsVars; j++) {
									LamAdd = (i == j) ? LamAbs : 0.0;
									buf[plainInd4Blk(i, j)] =
										0.5 * coef * (Jac_c_glob[i][j] + LamAdd);
								}
							}
							idxm = idGlob;
							idxn = idGlob;
							MatSetValuesBlocked(ctxKSP.A, 1, &idxm, 1, &idxn, buf, ADD_VALUES);
						}

						continue;
					}

					hsLogError(
						"t_DomEu1stImpl::makeTimeStep: unsupported face type for global Jac calculation");

				}	// iterate over cell faces

			}	// iterate over cells

		} // matrix & rhs insertions for the zone

	}	// for Zones

	VecAssemblyBegin(ctxKSP.b);
	VecAssemblyEnd(ctxKSP.b);

	MatAssemblyBegin(ctxKSP.A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ctxKSP.A, MAT_FINAL_ASSEMBLY);

	KSPSolve(ctxKSP.ksp, ctxKSP.b, ctxKSP.x);

	double norm;
	int nits;
	KSPGetResidualNorm(ctxKSP.ksp, &norm);
	KSPGetIterationNumber(ctxKSP.ksp, &nits);
	hsLogMessage("KSP: Norm of error %.6e iterations %d\n", norm, nits);

	// update my portion of solution
	double ResidTot = 0.0;
	double dU_max_Global = 0.0;
	{
		PetscScalar* ardU;  // processor's portion of the global field vector
		VecGetArray(ctxKSP.x, &ardU);

		int iRowBase;
		t_ConsVars dU;

		for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

			t_Zone& zne = Zones[iZone];

			for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

				iRowBase = getGlobInd(iZone, iCell) * NConsVars;
				for (int k = 0; k < NConsVars; k++)
					dU[k] = ardU[iRowBase + k];

				getCellCSV(iZone, iCell) += dU;

				double dU_norm = dU.norm();
				ResidTot += dU_norm;
				if (dU_norm > dU_max_Global) dU_max_Global = dU_norm;

			}
		}

		VecRestoreArray(ctxKSP.x, &ardU);
	}

	MPI_Allreduce(&ResidTot, &G_State.ResidTot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(&dU_max_Global, &dU_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (G_State.mpiRank == 0) {
		hsLogMessage("============");
		hsLogMessage("Time=%.6lf, Resid=%.6e, Local Max dU=%.6e", G_State.time, ResidTot, dU_max);
		hsLogMessage("Max normal velo resid at wall(sym):%.6e", G_State.ResidNormVeloWall);
	}

	G_State.time += dt;

	G_GhostMngEu.exchangeCSV();

	MatZeroEntries(ctxKSP.A);
	VecZeroEntries(ctxKSP.b);

}

void t_DomEu1stImpl::testKSP() {


	int idxm[3];
	int idxn[3];
	double vals[3];

	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

		t_Zone& zne = Zones[iZone];


		// insert some elements in matrix and probe ksp
		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			int idGlob = getGlobInd(iZone, iCell);
			int iRowBase = NConsVars * idGlob;

			for (int k = 0; k < NConsVars; k++) {
				if (k == 0) {
					idxm[0] = iRowBase + k;
					idxn[0] = iRowBase + k;
					idxn[1] = idxn[0] + 1;

					vals[0] = 1.001;
					vals[1] = -1;

					MatSetValues(ctxKSP.A, 1, idxm, 2, idxn, vals, ADD_VALUES);
					continue;
				}
				if (k == NConsVars - 1) {
					idxm[0] = iRowBase + k;
					idxn[0] = iRowBase + k - 1;
					idxn[1] = idxn[0] + 1;

					vals[0] = -1;
					vals[1] = 1.001;

					MatSetValues(ctxKSP.A, 1, idxm, 2, idxn, vals, ADD_VALUES);
					continue;
				}
				idxm[0] = iRowBase + k;
				idxn[0] = iRowBase + k - 1;
				idxn[1] = idxn[0] + 1;
				idxn[2] = idxn[1] + 1;
				vals[0] = -1;
				vals[1] = 2;
				vals[2] = -1;
				MatSetValues(ctxKSP.A, 1, idxm, 3, idxn, vals, ADD_VALUES);
			}
		}


	}

	MatAssemblyBegin(ctxKSP.A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ctxKSP.A, MAT_FINAL_ASSEMBLY);

	Vec exact_sol;
	VecCreate(PETSC_COMM_WORLD, &exact_sol);
	VecDuplicate(ctxKSP.x, &exact_sol);

	VecSet(exact_sol, 1.0);
	MatMult(ctxKSP.A, exact_sol, ctxKSP.b);

	KSPSolve(ctxKSP.ksp, ctxKSP.b, ctxKSP.x);

	VecAXPY(ctxKSP.x, -1.0, exact_sol);
	double norm;
	VecNorm(ctxKSP.x, NORM_2, &norm);
	int nits;
	KSPGetIterationNumber(ctxKSP.ksp, &nits);
	hsLogMessage("Norm of error %.6e iterations %d\n", norm, nits);

}