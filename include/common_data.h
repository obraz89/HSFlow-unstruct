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

static const int MaxNumVertsInFace = 4;

static const int MaxNumFacesInCell = 6;

static const int MaxNumVertsInCell = 8;

typedef long int lint;

template<int len>class t_Vec {
	double cont[len];
public:
	double& operator[](int i) {
#ifdef _DEBUG
		//if (i<0 || i>len-1) hsLogError("t_Vec:ind out of range");
#endif // _DEBUG
		return cont[i]; 
	}
	const double& operator[](int i) const { return cont[i]; }
};



template<int len> double dot(const t_Vec<len>& v1, const t_Vec<len>& v2) {
	double ret = 0.0;
	for (int i=0; i<len; i++) ret+=v1.cont[]
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

using t_Vec3d = t_Vec<3>;

enum struct t_FaceBCKind {
	Fluid = 0,
	Inflow,
	Outflow,
	Sym,
	Wall
};

enum struct t_CellKind {
	Tetra = 0,
	Brick

};

struct t_Cell;

struct t_Vert {

	lint Id;

	// list of cells that has this Vertex
	t_Cell *pNeibCells;

	t_Vec3d xyz;

};

struct t_Face {

	lint Id;

	int NumOfVerts;

	t_Vert Verts[MaxNumVertsInFace];

	t_Cell *pLeftCell, *pRightCell;

	// local indices of the face for left & right cells

	int IndLeftCellFace, IndRightCellFace;

	t_FaceBCKind BCKind;

	t_Vec3d Normal;

	t_Vec3d Center;

	void ComputeFaceCenter();

	void ComputeFaceNormal();

};



struct t_Cell {

	lint Id;

	int Nverts;

	int NFaces;

	t_Vert* Verts[MaxNumVertsInCell];
	t_Face* Faces[MaxNumFacesInCell];

	// +1 if Face Normal directed outward of the cell, -1 otherwise
	int FacesNormOutward[MaxNumFacesInCell];

	// List of Neighbor Cells
	int NumCellsNeig;
	t_Cell* CellsNeig;

	t_Vec3d Center;

	double Volume;

};

struct t_ZoneFacePatch {

	// Boundary condition on the face
	char szBC[33] = "";                 // BC-family name, MUST be empty if abutted

	// implement later
	//hsflow::TPhysBCCaps* BC = nullptr;  // NULL if abutted or not loaded yet

	bool isSkipped = false;     // face's grid layer skipped for processing by abutted zone

	t_Face* Faces;


};

struct t_Zone {

	char szName[40];  // name of the zone, initialized by '\0'

	lint nVerts;

	lint nCells;

	t_Vert *Verts;
	t_Cell *Cells;

};

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
