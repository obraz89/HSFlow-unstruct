#pragma once

#include "common_data.h"

static const int MaxNumVertsInFace = 4;

static const int MaxNumFacesInCell = 6;

static const int MaxNumVertsInCell = 8;

static const int MaxNumEdgesInCell = 12;

using t_BufFace2Vert = t_BufIndsStat<lint, MaxNumFacesInCell, MaxNumVertsInFace>;

using t_BufEdge2Vert = t_BufIndsStat<lint, MaxNumEdgesInCell, 2>;

// set of indexes, Face 2 Vertex
struct t_Vert;
using t_SetOfpVerts = t_TSet<const t_Vert*, MaxNumVertsInFace>;

struct t_Zone;

struct t_Cell;

struct t_Vert {

	lint Id;

	int NNeigCells;

	// list of cells that has this Vertex
	t_Cell** pNeigCells;

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
	t_Cell* pMyCell, * pOppCell;

	int IndLeftCellFace, IndRightCellFace;

	// 0 for fluid cell, others values are ids of bcs
	t_FaceBCID BCId;

	t_Vec3 Normal;

	t_Vec3 Center;

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

	int NCellsNeig() {
		int ret = 0;
		for (int i = 0; i < MaxNumFacesInCell; i++)
			if (pCellsNeig[i] != nullptr) ret++;
		return ret;
	};

	t_Cell() {
		for (int i = 0; i < MaxNumVertsInCell; i++) pVerts[i] = nullptr;
		for (int i = 0; i < MaxNumFacesInCell; i++) {
			pFaces[i] = nullptr;
			pCellsNeig[i] = nullptr;
		}
	}

	t_Vert& getVert(int ind) { return *pVerts[ind]; }
	const t_Vert& getVert(int ind) const { return *pVerts[ind]; }

	t_Vert* getpVert(int ind) { return pVerts[ind]; }
	const t_Vert* getpVert(int ind) const { return pVerts[ind]; }

	t_Face& getFace(int ind) { return *pFaces[ind]; }
	const t_Face& getFace(int ind) const { return *pFaces[ind]; }

	t_Face* getpFace(int ind) { return pFaces[ind]; }
	const t_Face* getpFace(int ind) const { return pFaces[ind]; }

	void setKind(t_CellKind a_Kind);
	void calcFaceNormalAreaOutward(int iface, t_Vec3& norm, double& area) const;

	bool isMyFace(int face_ind) const { return (this == pFaces[face_ind]->pMyCell); }

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


	// Global 0-based id of a zone on current worker
	int idGlob;
	// for now all vertices in the zone are real
	lint nVerts;
	// cells are stored as plain array
	lint nCellsReal;
	// total number of cells, real+ghosts
	lint nCellsTot;
	// number of faces in zone
	lint nFaces;
	// verts are real verts of the zone
	t_Vert* Verts;
	// plain array of cells, layout:
	// first nCellsReal are real cells
	// then comes ghostcells from zone0, then from zone1, etc...
	t_Cell* Cells;
	// real cells are decomposed into faces,
	// each face is stored once (no duplicates)
	t_Face* Faces;

public:

	void initialize(lint a_nVerts, lint a_nCellsReal, lint a_nCellsTot);

	void setIdGlob(int a_id) { idGlob = a_id; }
	int getIdGlob() const { return idGlob; }

	void setName(std::string a_name) { name = a_name; }
	const char* getName() const { return name.c_str(); }

	const lint& getnVerts() const { return nVerts; }
	const lint& getnCellsTot() const { return nCellsTot; }
	const lint& getnCellsReal() const { return nCellsReal; }
	lint getnCellsGhost() const { return nCellsTot - nCellsReal; }

	t_Cell& getCell(lint cell_ID) { return Cells[cell_ID]; }
	const t_Cell& getCell(lint cell_ID) const { return Cells[cell_ID]; }

	t_Cell* getpCell(lint cell_ID) { return &Cells[cell_ID]; }
	const t_Cell* getpCell(lint cell_ID) const { return &Cells[cell_ID]; }

	t_Vert& getVert(lint vert_ID) { return Verts[vert_ID]; }
	const t_Vert& getVert(lint vert_ID) const { return Verts[vert_ID]; }

	t_Vert* getpVert(lint vert_ID) { return &Verts[vert_ID]; }
	const t_Vert* getpVert(lint vert_ID) const { return &Verts[vert_ID]; }

	t_Face& getFace(lint face_ID) { return Faces[face_ID]; }
	const t_Face& getFace(lint face_ID) const { return Faces[face_ID]; }

	t_Face* getpFace(lint face_ID) { return &Faces[face_ID]; }
	const t_Face* getpFace(lint face_ID) const{ return &Faces[face_ID]; }

	void makeVertexConnectivity();
	void makeCellConnectivity();
	void makeFaces();
	void updateFacesWithBCPatch(const t_Face* face_patch, const int NFacesInPatch);

	std::vector<t_Cell*> getNeigCellsOfCellFace(const t_Cell& cell, int face_ind) const;
	void init_face2cell_conn(lint a_id, t_Cell& cell, int face_ind);

	int getFacePos(lint cell_id, const std::vector<lint> vert_ids) const;

	void getNeigAbutCellId(const std::vector<lint>& vert_ids, lint& cell_id, int& face_pos) const;

	bool isRealCell(lint cell_id) {
		// indices of real nodes are [0...nCellsReal-1]
		return (0 <= cell_id && cell_id <= nCellsReal - 1);
	}

	lint getNFaces() const { return nFaces; }

	~t_Zone() { delete[] Verts, Cells, Faces; }

};

struct t_Mesh
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
	//double gridCellScaleMin, gridCellScaleMax;
	//bool grid_update_distance_to_walls();
	//bool grid_make_symmetric_boco(const std::string& boco_name);

	void loadCells();

	void loadBCs();

	void makeVertexConnectivity();
	void makeCellConnectivity();
	void makeFaces();

	bool checkNormalOrientations();
	double calcUnitOstrogradResid();

	// Gas parameters
	//double(*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	//std::map<std::string, double> mapCasePrms_real;

	void initializeFromCtxStage1();

	void initializeFromCtxStage2();

	virtual void allocateFlowSolution() = 0;

	virtual void initializeFlow() = 0;

	// for debug
	virtual void dump_flow() = 0;
	virtual void dump_geom() = 0;

	virtual ~t_Mesh() {};
};

extern t_Mesh* G_pMesh;
