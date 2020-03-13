#pragma once

#include <map>

#include "dll_import-export.h"

//
// Problem solving state
//
struct TState
{
	int mpiRank, mpiNProcs;  // MPI rank, number of procs
							 // Relation of zone index and MPI rank working with it
	int* map_zone2rank;  // map_zone2rank[izne] == mpiRank, where izne -- 0-based zone index

						 //---

	int nTmStep;     // current time step number
	int nwtIter;     // current Newton's iteration number

					 /// L_inf norm of residual at 0-th Newton's iteration at current step
	double initialResidual;
};

DLLIMPEXP extern TState G_State;

struct t_Zone;

static const int MaxNumOfVertsInFace = 4;

static const int MaxNumberOfFacesInCell = 6;

typedef long int lint;

enum struct t_FaceBCKind {
	FluidFace = 0,
};

struct t_GrdVert {

	lint Id;

	double xyz[3];

};

struct t_GrdFace {

	lint Id;

	int NumOfVerts;

	t_GrdVert* Verts[MaxNumOfVertsInFace];

	lint IdLeftCell, IdRightCell;

	double Normal[3];

	double Center[3];

	void ComputeFaceCenter();

	void ComputeFaceNormal();

};



struct t_GrdCell {

	lint Id;

	t_GrdFace* Faces[MaxNumberOfFacesInCell];

};

struct t_Zone {};

struct DLLIMPEXP t_Domain
{
	int nu,	  // number of dependent (unknown) variables
		nDim; // number of independent variables (problem dimensions)

			  // Physical equations of the problem
			  //
	//hsflow::TPhysPluginBase* phys = nullptr;

	// Domain zones with grid and solution data
	//
	int nZones;  // total number of zones
	int bs, be;  // start & end (inclusive) 0-based zone (aka block) indices in the current MPI rank
	int* map_iZne2cgID;  // map_iZne2cgID[b] == cgZne, where b -- internal 0-based zone index, cgZne -- CGNS 1-based zone ID

	t_Zone* Zones;

	// Global grid info
	double gridCellScaleMin, gridCellScaleMax;
	bool grid_update_distance_to_walls();
	bool grid_make_symmetric_boco(const std::string& boco_name);

	// Gas parameters
	double(*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	std::map<std::string, double> mapCasePrms_real;
};

DLLIMPEXP extern t_Domain G_Domain;
