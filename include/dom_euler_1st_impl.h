#pragma once

#pragma once

#include "dom_euler_base.h"

#include "petsc.h"

// solving A*x = b
struct t_KSPSolver {

	Vec x, b;
	Mat A;
	KSP ksp;
	PetscReal resid;

	// size of vecs & mat
	int dim;

	// for debug only
	Vec exact_sol;

	~t_KSPSolver() {
		KSPDestroy(&ksp);
		VecDestroy(&x);
		VecDestroy(&b);  
		MatDestroy(&A);

		VecDestroy(&exact_sol);
	}


};

class t_DomEu1stImpl : public t_DomEuBase {
	t_KSPSolver ctxKSP;
public:
	void calcReconstData() {};
	void calcFaceFlux(int iZone, lint iFace);
	virtual void makeTimeStep();
	virtual void allocateFlowSolution();
	// local for this rank 0-based index of cell
	int getLocInd(const int iZone, const int iCell);
	// global (all-ranks) 0-based index of cell
	int getGlobInd(const int iZone, const int iCell);

	// test
	void makeTimeStep_SingleZone();
};

