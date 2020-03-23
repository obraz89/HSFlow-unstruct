#pragma once

#include <map>

#include "dll_import-export.h"

typedef __int64 lint;

static const int MaxNumVertsInFace = 4;

static const int MaxNumFacesInCell = 6;

static const int MaxNumVertsInCell = 8;

static const int MaxNumEdgesInCell = 12;

template<typename T>
class t_BufInds {
private:
	// TODO: for some reason, vs doesnt remove default constructor with
	// t_BufInds() = delete;
	// making it inaccessible old style
	t_BufInds() {};

	T* buf;
public:
	T nRows, nCols;
	T* data() { return buf; }
	//t_BufInds() = delete;
	t_BufInds(t_BufInds&) = delete;
	t_BufInds(T a_NR, T a_NC) { nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; }
	void allocate(T a_NR, T a_NC) { delete[] buf;  nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; };
	T& get_val(T i, T j) { return *(buf + i*nCols + j); };
	const T& get_val(T i, T j) const { return *(buf + i*nCols + j); }
	~t_BufInds() { delete[] buf; }

};
// container for small packs of index sets
template<typename T, int nRows, int nCols>
class t_BufIndsStat {
private:
	T buf[nRows][nCols];
public:
	T* data() { return buf; }
	T& get_val(int i, int j) { return buf[i][j]; };
	const T& get_val(int i, int j) const { return buf[i][j]; }
	~t_BufIndsStat() { delete[] buf; }
};

using t_BufFace2Vert = t_BufIndsStat<lint, MaxNumFacesInCell, MaxNumVertsInFace>;

using t_BufEdge2Vert = t_BufIndsStat<lint, MaxNumEdgesInCell, 2>;

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
	Brick,
	Pyra,
};

struct t_Cell;

struct t_Vert {

	lint Id;

	// list of cells that has this Vertex
	t_Cell *pNeigCells;

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

	t_CellKind Kind;

	int NVerts;

	int NEdges;

	int NFaces;

	t_Vert* pVerts[MaxNumVertsInCell];
	t_Face* pFaces[MaxNumFacesInCell];

	// +1 if Face Normal directed outward of the cell, -1 otherwise
	int FacesNormOutward[MaxNumFacesInCell];

	// List of Neighbor Cells
	int NumCellsNeig;
	t_Cell* CellsNeig;

	t_Vec3d Center;

	double Volume;

	t_Vert& getVert(int ind) { return *pVerts[ind]; }
	const t_Vert& getVert(int ind) const{ return *pVerts[ind]; }

	void setKind(t_CellKind a_Kind);

};

struct t_CellFaceList {

	const t_Cell& Cell;

	// array of faces via vertices
	t_BufFace2Vert F2V;

	// array containing number of vertexes for each face 
	int ArrNumOfVertsInFaces[MaxNumFacesInCell];

	int NFaces() const { return Cell.NFaces; };
	int NVertInFace(int indFace) const { return ArrNumOfVertsInFaces[indFace]; };

	t_CellFaceList() = delete;
	t_CellFaceList(const t_Cell& cell);

};

struct t_CellEdgeList {

	const t_Cell& Cell;

	// array of edges via vertices 
	t_BufEdge2Vert E2V;

	int NEdges() const { return Cell.NEdges; };

	t_CellEdgeList() = delete;
	t_CellEdgeList(const t_Cell& cell);



};

struct t_ZoneFacePatch {

	// Boundary condition on the face
	char szBC[33] = "";                 // BC-family name, MUST be empty if abutted

	// implement later
	//hsflow::TPhysBCCaps* BC = nullptr;  // NULL if abutted or not loaded yet

	bool isSkipped = false;     // face's grid layer skipped for processing by abutted zone

	t_Face* pFaces;


};

struct t_Zone {

	char szName[40];  // name of the zone, initialized by '\0'

	lint nVerts;

	lint nCells;

	t_Vert *pVerts;
	t_Cell *pCells;

	void initialize(lint nVerts, lint nCells);

	t_Cell& getCell(lint cell_ID) { return pCells[cell_ID]; }
	const t_Cell& getCell(lint cell_ID) const{ return pCells[cell_ID]; }

	t_Vert& getVert(lint vert_ID) { return pVerts[vert_ID]; }
	const t_Vert& getVert(lint vert_ID) const{ return pVerts[vert_ID]; }

	void getNeigCellsOfVertices();
	void makeFaceList();

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

	void makeFaceLists();

	// Gas parameters
	double(*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	std::map<std::string, double> mapCasePrms_real;
};

DLLIMPEXP extern t_Domain G_Domain;
