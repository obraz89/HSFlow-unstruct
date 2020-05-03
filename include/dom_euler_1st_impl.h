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

	// number of rows that are stored at this rank
	int dimMy;
	// number of rows of global matrix
	int dimGlob;

	~t_KSPSolver() {
		KSPDestroy(&ksp);
		VecDestroy(&x);
		VecDestroy(&b);  
		MatDestroy(&A);

	}


};

// absolute value of largest eigenvalue of the linearized Riemann problem between cells c and d
// i.e. store this for every face
struct t_ZoneLambdasCD {
	double* Lambdas;
	~t_ZoneLambdasCD() { delete[] Lambdas; }
};

class t_DomEu1stImpl : public t_DomEuBase {
	t_ZoneLambdasCD* ZonesLambdaCD;
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
	void testKSP();
	virtual ~t_DomEu1stImpl() {
		for (int i = iZneMPIs; i <= iZneMPIe; i++) {
			delete[] ZonesLambdaCD[i].Lambdas;
		}
		delete[] ZonesLambdaCD;
	}
};

