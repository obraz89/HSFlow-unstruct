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

	// get dimensions of global unknown vector du
	// and starting offset for each rank
	int NCellsMy = 0;

	for (int i = iZneMPIs; i <= iZneMPIe; i++)
			NCellsMy += Zones[i].getnCellsReal();

	ctxKSP.dim = NCellsMy*NConsVars;

	// number of non-zeros in diagonal (d_nnz) and non-diagonal (o_nnz) 
	// part for each row in ownership of this rank
	int* d_nnz = new int[ctxKSP.dim];
	int* o_nnz = new int[ctxKSP.dim];

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
	// allocate vectors
	perr = VecCreate(PETSC_COMM_WORLD, &ctxKSP.x);
	perr = VecSetSizes(ctxKSP.x, ctxKSP.dim, PETSC_DECIDE);
	perr = VecSetBlockSize(ctxKSP.x, NConsVars);
	perr = VecSetType(ctxKSP.x, VECMPI);
	if (perr) hsTHROW("Failed to init PETSc field vector");

	VecDuplicate(ctxKSP.x, &ctxKSP.b);
	VecDuplicate(ctxKSP.x, &ctxKSP.exact_sol);

	// allocate matrix
	int NRowsMy = NCellsMy * NConsVars;

	perr = MatCreateAIJ(PETSC_COMM_WORLD,
		NRowsMy, NRowsMy,
		PETSC_DETERMINE, PETSC_DETERMINE,
		0, d_nnz, 0, o_nnz, &ctxKSP.A);
	if (perr) hsTHROW("Failed to init PETSc Jacobi matrix");

	//MatSetOption(ctx.matJac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

	perr = KSPCreate(PETSC_COMM_WORLD, &ctxKSP.ksp);
	// deflated GMRES, the adaptive strategy allows to switch to the deflated GMRES when the stagnation occurs
	//perr = KSPSetType(ctxKSP.ksp, KSPDGMRES); 

	KSPSetOperators(ctxKSP.ksp, ctxKSP.A, ctxKSP.A);
	// ???
	KSPSetTolerances(ctxKSP.ksp, 1.e-5, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);

	KSPSetFromOptions(ctxKSP.ksp);
	
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



}

void t_DomEu1stImpl::makeTimeStep() {

	hsLogMessage("Hello world");

	double dt_local = calcDt();

	double dt;

	MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	double dU_max = 0.0;

	hsLogMessage("Computed dt=%lf", dt);

	G_State.ResidTot = 0.0;

	G_State.ResidNormVeloWall = 0.0;

	int idxm;
	int idxn;
	double val = 1.0;
	for (int iZone = iZneMPIs; iZone <= iZneMPIe; iZone++) {

		t_Zone& zne = Zones[iZone];

		//for (int iFace = 0; iFace < zne.getNFaces(); iFace++) {
		//	calcFaceFlux(iZone, iFace);
		//}

		// testing : inserting just diag elems
		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			int idGlob = getGlobInd(iZone, iCell);
			int iRowBase = NConsVars * idGlob;

			for (int k = 0; k < NConsVars; k++) {
				idxm = iRowBase + k;
				idxn = iRowBase + k;
				MatSetValues(ctxKSP.A, 1, &idxm, 1, &idxn, &val, ADD_VALUES);
			}
		}


	}

	MatAssemblyBegin(ctxKSP.A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ctxKSP.A, MAT_FINAL_ASSEMBLY);

	VecSet(ctxKSP.exact_sol, 1.0);
	MatMult(ctxKSP.A, ctxKSP.exact_sol, ctxKSP.b);

	KSPSolve(ctxKSP.ksp, ctxKSP.b, ctxKSP.x);

	VecAXPY(ctxKSP.x, -1.0, ctxKSP.exact_sol);
	double norm;
	VecNorm(ctxKSP.x, NORM_2, &norm);
	int nits;
	KSPGetIterationNumber(ctxKSP.ksp, &nits);

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