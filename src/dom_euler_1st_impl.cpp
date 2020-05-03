#include "dom_euler_1st_impl.h"

#include "flux_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "ghost_euler.h"

int t_DomEu1stImpl::getLocInd(const int iZone, const int iCell) {

	int idZoneOffset = 0;

	for (int i = iZneMPIs; i < iZone; i++)
		idZoneOffset += Zones[i].getnCellsReal();

	return idZoneOffset + iCell;


}

int t_DomEu1stImpl::getGlobInd(const int iZone, const int iCell) {

	int idZoneOffset = 0;

	for (int i = 0; i < iZone; i++)
		idZoneOffset += Zones[i].getnCellsReal();

	return idZoneOffset + iCell;

};

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

	// number of non-zeros in diagonal (d_nnz) and non-diagonal (o_nnz) 
	// part for each row in ownership of this rank
	int* d_nnz = new int[ctxKSP.dimMy];
	int* o_nnz = new int[ctxKSP.dimMy];

	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {
		const t_Zone& Zne = Zones[iZone];
		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {
			int idMy = getLocInd(iZone, iCell);
			// block non-zeros are the same
			int dz = NConsVars;
			int oz = Zne.getCell(iCell).NCellsNeig() * NConsVars;
			for (int k = 0; k < NConsVars; k++) {
				int iRow = NConsVars * idMy + k;
				d_nnz[iRow] = dz;
				o_nnz[iRow] = oz;
			}
		}
	}

	PetscErrorCode perr;
	// create vectors
	// x - let petsc decide storage
	perr = VecCreate(PETSC_COMM_WORLD, &ctxKSP.x);
	perr = VecSetSizes(ctxKSP.x, PETSC_DECIDE, ctxKSP.dimGlob);
	if (perr) hsLogMessage("Error:Failed to init PETSc solution vector");

	// create rhs
	perr = VecCreateMPI(PETSC_COMM_WORLD, ctxKSP.dimMy, PETSC_DETERMINE, &ctxKSP.b);
	if (perr) hsLogMessage("Error:Failed to init PETSc rhs vector");


	// allocate matrix

	perr = MatCreateAIJ(PETSC_COMM_WORLD,
		ctxKSP.dimMy, ctxKSP.dimMy,
		PETSC_DETERMINE, PETSC_DETERMINE,
		0, d_nnz, 0, o_nnz, &ctxKSP.A);
	if (perr) hsLogMessage("Error:Failed to init PETSc matrix");

	//MatSetOption(ctx.matJac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

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

void t_DomEu1stImpl::makeTimeStep_SingleZone() {

	Vec v;

}

void t_DomEu1stImpl::calcFaceFlux(int iZone, lint iFace) {

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

					const t_Flux& flux = getFlux(iZone, face.Id);

					dU_rhs += coef * flux;

				}

				// insert values in my portion of rhs vector
				for (int k = 0; k < NConsVars; k++) {
					buf[k] = dU_rhs[k];
					idx[k] = iRowBase + k;
				}

				VecSetValues(ctxKSP.b, NConsVars, idx, buf, INSERT_VALUES);

			}

			VecAssemblyBegin(ctxKSP.b);
			VecAssemblyEnd(ctxKSP.b);
		
		}



		// constructing The Matrix 
		{
			for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

				t_Cell& cell = zne.getCell(iCell);

				for (int j = 0; j < cell.NFaces; j++) {

				}

			}

			MatAssemblyBegin(ctxKSP.A, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(ctxKSP.A, MAT_FINAL_ASSEMBLY);

		}

	}

	hsLogMessage("dfgsdfgsdfg");
	return;

	KSPSolve(ctxKSP.ksp, ctxKSP.b, ctxKSP.x);

	double norm;
	int nits;
	KSPGetResidualNorm(ctxKSP.ksp, &norm);
	KSPGetIterationNumber(ctxKSP.ksp, &nits);
	hsLogMessage("KSP: Norm of error %.6e iterations %d\n", norm, nits);

	return;


	double ResidTot;
	MPI_Allreduce(&G_State.ResidTot, &ResidTot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double dU_max_Global;
	MPI_Allreduce(&dU_max, &dU_max_Global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (G_State.mpiRank == 0) {
		hsLogMessage("============");
		hsLogMessage("Time=%.6lf, Resid=%.6e, Local Max dU=%.6e", G_State.time, ResidTot, dU_max);
		hsLogMessage("Max normal velo resid at wall(sym):%.6e", G_State.ResidNormVeloWall);
	}

	G_State.time += dt;

	G_GhostMngEu.exchangeCSV();

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