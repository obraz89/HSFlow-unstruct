#pragma once

#include <map>

#include "dll_import-export.h"

#include <array>
#include <vector>

#include "logging.h"

#include "matrix_small.h"

#include "flow_model.h"

typedef __int64 lint;

static const int MaxNumVertsInFace = 4;

static const int MaxNumFacesInCell = 6;

static const int MaxNumVertsInCell = 8;

static const int MaxNumEdgesInCell = 12;

template<typename T>
class t_BufInds {
	T* buf;
public:
	T nRows, nCols;
	T* data() { return buf; }
	const T* data() const { return buf; }
	t_BufInds() = delete;
	t_BufInds(t_BufInds&) = delete;
	t_BufInds(T a_NR, T a_NC) { nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; }
	void allocate(T a_NR, T a_NC) { delete[] buf;  nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; };
	T& get_val(T i, T j) { return *(buf + i*nCols + j); };
	const T& get_val(T i, T j) const { return *(buf + i*nCols + j); }
	~t_BufInds() { delete[] buf; }

};

// plain set of indices, usually decomposition of a face
template<typename T, int Nmax>
class t_TSet {
	std::array<T, Nmax> buf;
	int NElems;
public:
	t_TSet() :NElems(0) { for (int i = 0; i < Nmax; i++) buf[i] = 0; }
	int size() const { return NElems; }
	void setSize(int a_size) { 
#ifdef _DEBUG
		if (a_size > Nmax) hsLogMessage("t_TSet:setSize:Error: size is too big");
#endif // _DEBUG
		NElems = a_size;
	}
	T& operator[](int i) { 
#ifdef _DEBUG
		if ((i > Nmax-1) || (i<0)) hsLogMessage("t_TSet: wrong index");
#endif // _DEBUG
		return buf[i]; }
	const T& operator[](int i) const { 
#ifdef _DEBUG
		if ((i > Nmax - 1) || (i<0)) hsLogMessage("t_TSet: wrong index");
#endif // _DEBUG
		return buf[i]; }

	// weak comparison, order not important {1,2,3,4}=={1,3,4,2} : true
	static bool cmp_weak(const t_TSet<T, Nmax>& lv, const t_TSet<T, Nmax>& rv){
		if (lv.size() != rv.size()) return false;
		std::array<T, Nmax> l_sorted = lv.buf;
		std::sort(l_sorted.begin(), l_sorted.end());
		std::array<T, Nmax> r_sorted = rv.buf;
		std::sort(r_sorted.begin(), r_sorted.end());
		bool cmp_ok = true;
		for (int i = 0; i < Nmax; i++) cmp_ok = cmp_ok && (l_sorted[i] == r_sorted[i]);
		return cmp_ok;
	};

	static bool cmp_strict(const t_TSet<T, Nmax>& lv, const t_TSet<T, Nmax>& rv){
		if (lv.size() != rv.size()) return false;
		bool cmp_ok = true;
		for (int i = 0; i < NMax; i++) cmp_ok = cmp_ok && (l_sorted[i] == r_sorted[i]);
		return cmp_ok;
	}
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
};

using t_BufFace2Vert = t_BufIndsStat<lint, MaxNumFacesInCell, MaxNumVertsInFace>;

using t_BufEdge2Vert = t_BufIndsStat<lint, MaxNumEdgesInCell, 2>;

// set of indexes, Face 2 Vertex
struct t_Vert;
using t_SetOfpVerts = t_TSet<const t_Vert*, MaxNumVertsInFace>;

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

	int NNeigCells;

	// list of cells that has this Vertex
	t_Cell **pNeigCells;

	t_Vec3 xyz;

	~t_Vert() { delete[] pNeigCells; }

};

struct t_Face {

	lint Id;

	int NVerts;

	const t_Vert* pVerts[MaxNumVertsInFace];

	// MyCell: Face coincide with cell face, scalar product of normals +1
	// OppCell: opposing cell, scalar product of their normals -1 
	// MyCell is always defined but OppCell can be missing (BC face)
	t_Cell *pMyCell, *pOppCell;

	int IndLeftCellFace, IndRightCellFace;

	t_FaceBCKind BCKind;

	t_Vec3 Normal;

	t_Vec3 Center;

	t_Flux Flux;

	double Area;

	t_Face() { pMyCell = nullptr; pOppCell = nullptr; }

	void ComputeFaceCenter();

	void ComputeFaceNormal();

	t_SetOfpVerts getVerts() const;

};



struct t_Cell {

	lint Id;

	t_CellKind Kind;

	int NVerts;

	int NEdges;

	int NFaces;

	t_Vert* pVerts[MaxNumVertsInCell];
	t_Face* pFaces[MaxNumFacesInCell];

	t_Cell* pCellsNeig[MaxNumFacesInCell];
	// store index of face for neighbor cell that corresponds to the particular face
	// this is to reduce computations
	int FaceIndNeig[MaxNumFacesInCell];
	// Number of Neighbor Fluid Cells

	t_Vec3 Center;

	double Volume;

	// Flow Data

	t_PrimVars PrimVars;

	int NCellsNeig() { int ret=0; 
		for (int i = 0; i < MaxNumFacesInCell; i++) 
			if (pCellsNeig[i] != nullptr) ret++; 
		return ret; 
	};

	t_Cell() {
		for (int i = 0; i < MaxNumVertsInCell;i++) pVerts[i] = nullptr;
		for (int i = 0; i < MaxNumFacesInCell; i++) {
			pFaces[i] = nullptr;
			pCellsNeig[i] = nullptr;
		}
	}

	t_Vert& getVert(int ind) { return *pVerts[ind]; }
	const t_Vert& getVert(int ind) const{ return *pVerts[ind]; }

	t_Vert* getpVert(int ind) { return pVerts[ind]; }
	const t_Vert* getpVert(int ind) const { return pVerts[ind]; }

	t_Face& getFace(int ind) { return *pFaces[ind]; }
	const t_Face& getFace(int ind) const{ return *pFaces[ind]; }

	t_Face* getpFace(int ind) { return pFaces[ind]; }
	const t_Face* getpFace(int ind) const{ return pFaces[ind]; }

	void setKind(t_CellKind a_Kind);
	void calcFaceNormalAreaOutward(int iface, t_Vec3& norm, double& area) const;

	bool isMyFace(int face_ind) const{ return (this == pFaces[face_ind]->pMyCell); }
	
	t_Vec3 getFaceNormalOutward(int iface) const;

};

struct t_CellFaceList {

	const t_Cell* pCell;

	t_SetOfpVerts F2V[MaxNumFacesInCell];
	// array of faces via vertices
	
	int NFaces() const { return pCell->NFaces; };
	int NVertInFace(int indFace) const { return F2V[indFace].size(); };

	t_CellFaceList() = delete;
	t_CellFaceList(const t_Cell& cell);

	const t_SetOfpVerts& getVertices(int indFace) const;

};

struct t_CellEdgeList {

	const t_Cell* pCell;

	// array of edges via vertices 
	t_BufEdge2Vert E2V;

	int NEdges() const { return pCell->NEdges; };

	t_CellEdgeList() :pCell(nullptr) {};
	void init(const t_Cell& cell);



};

struct t_ZoneFacePatch {

	// Boundary condition on the face
	char szBC[33] = "";                 // BC-family name, MUST be empty if abutted

	// implement later
	//hsflow::TPhysBCCaps* BC = nullptr;  // NULL if abutted or not loaded yet

	bool isSkipped = false;     // face's grid layer skipped for processing by abutted zone

	t_Face* pFaces;


};

class t_Zone {

	std::string name;  // name of the zone, initialized by '\0'

	int id;

	lint nVerts;

	lint nCells;

	lint nFaces;

	t_Vert *Verts;
	t_Cell *Cells;
	t_Face *Faces;

public:

	void initialize(lint nVerts, lint nCells);

	void setId(int a_id) { id = a_id; }

	void setName(const char* a_name) { name = std::string(a_name); }
	const char* getName() const { return name.c_str(); }

	const lint& getnVerts() const { return nVerts; }
	const lint& getnCells() const { return nCells; }

	t_Cell& getCell(lint cell_ID) { return Cells[cell_ID]; }
	const t_Cell& getCell(lint cell_ID) const{ return Cells[cell_ID]; }

	t_Cell* getpCell(lint cell_ID) { return &Cells[cell_ID]; }
	const t_Cell* getpCell(lint cell_ID) const{ return &Cells[cell_ID]; }

	t_Vert& getVert(lint vert_ID) { return Verts[vert_ID]; }
	const t_Vert& getVert(lint vert_ID) const{ return Verts[vert_ID]; }

	t_Vert* getpVert(lint vert_ID) { return &Verts[vert_ID]; }
	const t_Vert* getpVert(lint vert_ID) const{ return &Verts[vert_ID]; }

	void makeVertexConnectivity();
	void makeCellConnectivity();
	void makeFaces();
	void updateFacesWithBCPatch(const t_Face* face_patch, const int NFacesInPatch);

	std::vector<t_Cell*> getNeigCellsOfCellFace(const t_Cell& cell, int face_ind) const;
	void init_face2cell_conn(lint a_id, t_Cell& cell, int face_ind);

	// flow solver

	double calcDt();
	void makeTimeStep(double dt);
	void calcFaceFlux(lint iFace);

	~t_Zone() { delete[] Verts, Cells, Faces; }

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

	void makeVertexConnectivity();
	void makeCellConnectivity();
	void makeFaces();

	bool checkNormalOrientations();
	double calcUnitOstrogradResid();

	void makeTimeStep();

	// Gas parameters
	double(*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	std::map<std::string, double> mapCasePrms_real;
};

DLLIMPEXP extern t_Domain G_Domain;
