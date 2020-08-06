#pragma once

#include "common_data_unstruct.h"

static const int MaxNumVertsInFace = 4;

static const int MaxNumFacesInCell = 6;

static const int MaxNumVertsInCell = 8;

static const int MaxNumEdgesInCell = 12;

// max number of tetras when making decomposition of a cell
static const int MaxNumTetrasFromCell = 12;

using t_BufFace2Vert = t_BufIndsStat<lint, MaxNumFacesInCell, MaxNumVertsInFace>;

using t_BufEdge2Vert = t_BufIndsStat<lint, MaxNumEdgesInCell, 2>;

// set of indexes, Face 2 Vertex
struct t_Vert;
using t_SetOfpVerts = t_TSet<const t_Vert*, MaxNumVertsInFace>;

class t_Zone;

struct t_Cell;

struct t_Vert {

	lint Id;

	int NNeigCells;

	// list of cells that has this Vertex
	t_Cell** pNeigCells;

	// weights of adjacent cells to get vertex flow vars
	double* pNeigCoefs;

	void calcAllocNeigCoefs();

	t_Vec3 xyz;
	t_Vert():pNeigCells(nullptr), pNeigCoefs(nullptr), Id(-1), NNeigCells(0) {}
	~t_Vert() { 
		if (pNeigCells !=nullptr) delete[] pNeigCells; 
		if (pNeigCoefs != nullptr) delete[] pNeigCoefs;
	}

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

	// index of vertex with max angle 
	// it is the root to compute base vectors in plane
	int IndVertRoot;

	// matrix for gradient reconstruction of flow vars
	t_SqMat3 MatGrad;

	// for debug
	double AngEdgeMin, AngEdgeMax;

	t_Face() { pMyCell = nullptr; pOppCell = nullptr; }

	void ComputeFaceCenter();

	void ComputeRootVertex();

	void ComputeMatGrad();

	bool isFluid() const { return BCId.get() == t_FaceBCID::Fluid; }

	void _makeCyclicListofVerts(const t_Vert* (&lst)[MaxNumVertsInFace + 2]) const;

	// reconstruct face gradients for arbitrary set of values
	// for example, in NS computations it is convinient to calculate for RUVWT
	// IMPORTANT: finite differences must be consistent with gradient matrix!
	template<int NVars> void ComputeFaceGrad(const t_Vec<NVars>& Umy, const t_Vec<NVars>& Uop, 
		const t_Mat<NVars, MaxNumVertsInFace>& UVerts, t_Mat<3, NVars>& Grads) {

		// make vertex values cyclic
		t_Mat<NVars, MaxNumVertsInFace + 2> UVertsCyc;

		for (int ivar = 0; ivar < NVars; ivar++) {

			UVertsCyc[ivar][0] = UVerts[ivar][NVerts - 1];

			for (int j = 0; j < NVerts; j++) {

				int jcyc = j + 1;

				UVertsCyc[ivar][jcyc] = UVerts[ivar][j];

			}

			UVertsCyc[ivar][NVerts + 1] = UVerts[ivar][0];

		}

		t_Vec3 RHS, grad;

		for (int ivar = 0; ivar < NVars; ivar++) {

			RHS[0] = Uop[ivar] - Umy[ivar];

			if (NVerts == 3) {

				int jcyc = IndVertRoot + 1;

				RHS[1] = UVertsCyc[ivar][jcyc + 1] - UVertsCyc[ivar][jcyc];
				RHS[2] = UVertsCyc[ivar][jcyc - 1] - UVertsCyc[ivar][jcyc];
			}

			if (NVerts == 4) {

				RHS[1] = 0.5*(UVerts[ivar][2] + UVerts[ivar][3] -
					          UVerts[ivar][0] - UVerts[ivar][1]);

				RHS[2] = 0.5*(UVerts[ivar][1] + UVerts[ivar][2] -
					          UVerts[ivar][3] - UVerts[ivar][0]);

			}

			grad = MatGrad*RHS;

			Grads.setCol(ivar, grad);

		}

	};

	t_SetOfpVerts getVerts() const;

};



struct t_Cell {

	lint Id;

	t_CellKind Kind;

	int NVerts;

	int NEdges;

	int NFaces;

	const t_Vert* pVerts[MaxNumVertsInCell];
	const t_Face* pFaces[MaxNumFacesInCell];

	t_Cell* pCellsNeig[MaxNumFacesInCell];
	// store index of face for neighbor cell that corresponds to the particular face
	// this is to reduce computations
	int FaceIndNeig[MaxNumFacesInCell];
	// Number of Neighbor Fluid Cells

	t_Vec3 Center;

	double Volume;

	double Diameter;

	int NCellsNeig() const{
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

	//t_Vert& getVert(int ind) { return *pVerts[ind]; }
	const t_Vert& getVert(int ind) const { return *pVerts[ind]; }

	//t_Vert* getpVert(int ind) { return pVerts[ind]; }
	const t_Vert* getpVert(int ind) const { return pVerts[ind]; }

	//t_Face& getFace(int ind) { return *pFaces[ind]; }
	const t_Face& getFace(int ind) const { return *pFaces[ind]; }

	//t_Face* getpFace(int ind) { return pFaces[ind]; }
	const t_Face* getpFace(int ind) const { return pFaces[ind]; }

	void setKind(t_CellKind a_Kind);
	void calcFaceNormalAreaOutward(int iface, t_Vec3& norm, double& area) const;

	bool isMyFace(int face_ind) const { return (this == pFaces[face_ind]->pMyCell); }

	t_Vec3 getFaceNormalOutward(int iface) const;

	void computeCenter();
	void computeDiameter();
	void computeVolume();

};

// decomposition of cell into tetras
class t_CellTetraList {

	t_Cell Tetras[MaxNumTetrasFromCell];
	// tmp verts neede for calcs, like cell center
	std::vector<t_Vert> VertsTmp;
	int Size;
public:
	int getSize() const { return Size; }
	t_Cell& getTetra(int ind) { return Tetras[ind]; }
	t_CellTetraList(const t_Cell& cell);

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
	t_SetOfpVerts& getVertices(int indFace);

};

struct t_CellEdgeList {

	const t_Cell* pCell;

	// array of edges via vertices 
	t_BufEdge2Vert E2V;

	int NEdges() const { return pCell->NEdges; };

	t_CellEdgeList() :pCell(nullptr) {};
	t_CellEdgeList(const t_Cell& cell) { init(cell); }
	void init(const t_Cell& cell);



};
// cells are stored as continuous range of cells with same kind
// layout must be the same as in original mesh
struct t_CellKindRange {
	t_CellKind kind;
	lint idStart;
	lint idEnd;
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
	int idGlob = -1;
	// for now all vertices in the zone are real
	lint nVerts = 0;
	// cells are stored as plain array
	lint nCellsReal = 0;
	// total number of cells, real+ghosts
	lint nCellsTot = 0;
	// number of faces in zone
	lint nFaces = 0;
	// number of BC faces
	lint nFacesBC = 0;
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
	// for zones we do not own, keep sizes
	void setNVertsNCells(lint a_nVerts, lint a_nCellsReal) {
		nVerts = a_nVerts;
		nCellsReal = a_nCellsReal;
	}

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

	std::vector<t_CellKindRange> getCellsOffsets() const;

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

	void makeVertexConnectivityRealCells();

	void makeCellConnectivity();
	void makeFaces();
	void updateFacesWithBCPatch(const t_Face* face_patch, const int NFacesInPatch);

	std::vector<t_Cell*> getNeigCellsOfCellFace(const t_Cell& cell, int face_ind) const;
	void init_face2cell_conn(lint a_id, t_Cell& cell, int face_ind);

	int getFacePos(lint cell_id, const std::vector<lint> vert_ids) const;

	void getNeigAbutCellId(const std::vector<lint>& vert_ids, lint& cell_id, int& face_pos) const;

	bool isRealCell(lint cell_id) const{
		// indices of real nodes are [0...nCellsReal-1]
		return (0 <= cell_id && cell_id <= nCellsReal - 1);
	}

	lint getNFaces() const { return nFaces; }
	lint getNFacesBC() const { return nFacesBC; }

	~t_Zone() { delete[] Verts, Cells, Faces; }

};

struct t_DomBase
{

	const int nDim = 3; // number of independent variables (problem dimensions)

	// Physical equations of the problem
	//
	//hsflow::TPhysPluginBase* phys = nullptr;

	// Domain zones with grid and solution data
	//
	int nZones;  // total number of zones
	int iZneMPIs, iZneMPIe;  // start & end (inclusive) 0-based zone (aka block) indices in the current MPI rank
	int* map_iZne2cgID;  // map_iZne2cgID[b] == cgZne, where b -- internal 0-based zone index, cgZne -- CGNS 1-based zone ID

	t_Zone* Zones;

	int nZonesMy() const { return iZneMPIe - iZneMPIs + 1; }

	// Global grid info
	//double gridCellScaleMin, gridCellScaleMax;
	//bool grid_update_distance_to_walls();
	//bool grid_make_symmetric_boco(const std::string& boco_name);

	void loadCells();

	void makeVertexConnectivityRealCells();

	void makeCellConnectivity();
	void makeFaces();

	bool checkNormalOrientations();
	double calcUnitOstrogradResid();

	// Gas parameters
	//double(*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	//std::map<std::string, double> mapCasePrms_real;

	void initializeFromCtxStage1();
	bool assignZonesToProcs();

	void initializeFromCtxStage2();

	virtual void checkMesh();

	void calcCellWeightsForVertices();

	// interface for the flow domain
	virtual void allocateFlowSolution() = 0;
	virtual void initializeFlow() = 0;
	virtual void prepareBeforeTimeMarch() = 0;
	virtual double loadField(std::string FieldName) = 0;
	virtual int getNu() const= 0;
	virtual std::vector<std::string> getFuncNamesIO() const = 0;
	virtual void getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const= 0;
	virtual void checkMinMaxCSV() = 0;
	virtual void makeTimeStep() = 0;

	virtual ~t_DomBase() { delete[] Zones; };
};

extern t_DomBase* g_pDomUnst;
