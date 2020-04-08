#include "mesh.h"

// TODO: put ghost manager inside mesh ?
#include "ghost_manager.h"

#include <fstream>

#include "CGNS-ctx.h"

// TODO: remove dependence from model by inheritances
#include "bc_data.h"

t_Mesh* G_pMesh;

t_CellKind getElementKind(CG_ElementType_t cg_type) {

	t_CellKind cell_kind = t_CellKind::None;

	if (cg_type == CG_TETRA_4) {
		cell_kind = t_CellKind::Tetra;
	}
	if (cg_type == CG_HEXA_8) {
		cell_kind = t_CellKind::Brick;
	}
	return cell_kind;

};

//************************************* Cell methods
void t_Cell::setKind(t_CellKind a_Kind) {

	Kind = a_Kind;

	if (Kind == t_CellKind::Brick) {

		NVerts = 8;

		NEdges = 12;

		NFaces = 6;

		return;
	}

	if (Kind == t_CellKind::Tetra) {

		NVerts = 4;

		NEdges = 6;

		NFaces = 4;

		return;
	}

	hsLogMessage("t_Cell::setSizes: can't handle this element yet");

};

// If cell vertices are read according to cgns vertex numbering,
// the computed normal is directed outward of the cell
// TODO: if grids with non-standard vertex numbering are read
// additional work should be done to ensure the normal is directed outward of the cell
void t_Cell::calcFaceNormalAreaOutward(int iface, t_Vec3& norm, double& area) const {

	t_CellFaceList flist(*this);

	// compute directly cell face normal
	t_SetOfpVerts verts = flist.getVertices(iface);

	if (verts.size() == 3) {
		t_Vec3 pnts[3];
		for (int k = 0; k < 3; k++) pnts[k] = verts[k]->xyz;
		ComputeTriangleAreaNormal(pnts, norm, area);
	}
	if (verts.size() == 4) {
		t_Vec3 pnts[4];
		for (int k = 0; k < 4; k++) pnts[k] = verts[k]->xyz;
		ComputeQuadAreaNormal(pnts, norm, area);

	}

	return;

}

t_Vec3 t_Cell::getFaceNormalOutward(int iface) const {

	t_Vec3 norm;

	const t_Face* pface = pFaces[iface];

	if (this == pface->pMyCell) norm = pface->Normal;
	else norm = -1.0 * pface->Normal;

	return norm;
}

void t_CellEdgeList::init(const t_Cell& a_Cell) {

	pCell = &a_Cell;

	if (pCell->Kind == t_CellKind::Brick) {

		// list of edges for a hexa cell
		// decomposition according cgns documentation: sids/conv.html#unst_hexa
		// vertexes (1,2,3,4,5,6,7,8) => edges 
		// (1,2), (2,3), (3,4), (4,1),
		// (1,5), (2,6), (3,7), (4,8),
		// (5,6), (6,7), (7,8), (8,5)

		const lint& V1 = pCell->getVert(0).Id;
		const lint& V2 = pCell->getVert(1).Id;
		const lint& V3 = pCell->getVert(2).Id;
		const lint& V4 = pCell->getVert(3).Id;
		const lint& V5 = pCell->getVert(4).Id;
		const lint& V6 = pCell->getVert(5).Id;
		const lint& V7 = pCell->getVert(6).Id;
		const lint& V8 = pCell->getVert(7).Id;

		E2V.get_val(0, 0) = V1;
		E2V.get_val(0, 1) = V2;

		E2V.get_val(1, 0) = V2;
		E2V.get_val(1, 1) = V3;

		E2V.get_val(2, 0) = V3;
		E2V.get_val(2, 1) = V4;

		E2V.get_val(3, 0) = V4;
		E2V.get_val(3, 1) = V1;

		E2V.get_val(4, 0) = V1;
		E2V.get_val(4, 1) = V5;

		E2V.get_val(5, 0) = V2;
		E2V.get_val(5, 1) = V6;

		E2V.get_val(6, 0) = V3;
		E2V.get_val(6, 1) = V7;

		E2V.get_val(7, 0) = V4;
		E2V.get_val(7, 1) = V8;

		E2V.get_val(8, 0) = V5;
		E2V.get_val(8, 1) = V6;

		E2V.get_val(9, 0) = V6;
		E2V.get_val(9, 1) = V7;

		E2V.get_val(10, 0) = V7;
		E2V.get_val(10, 1) = V8;

		E2V.get_val(11, 0) = V8;
		E2V.get_val(11, 1) = V5;

		return;
	}

	if (pCell->Kind == t_CellKind::Tetra) {

		// list of edges for a tetra cell
		// decomposition according cgns documentation: sids/conv.html#unst_tetra
		// vertexes (1,2,3,4) => edges (1,2),(2,3),(3,1),(1,4),(2,4),(3,4)

		const lint& V1 = pCell->getVert(0).Id;
		const lint& V2 = pCell->getVert(1).Id;
		const lint& V3 = pCell->getVert(2).Id;
		const lint& V4 = pCell->getVert(3).Id;

		E2V.get_val(0, 0) = V1;
		E2V.get_val(0, 1) = V2;

		E2V.get_val(1, 0) = V2;
		E2V.get_val(1, 1) = V3;

		E2V.get_val(2, 0) = V3;
		E2V.get_val(2, 1) = V1;

		E2V.get_val(3, 0) = V1;
		E2V.get_val(3, 1) = V4;

		E2V.get_val(4, 0) = V2;
		E2V.get_val(4, 1) = V4;

		E2V.get_val(5, 0) = V3;
		E2V.get_val(5, 1) = V4;

		return;
	}
	hsLogMessage("t_CellEdgeList: unsupported element type");
};

t_CellFaceList::t_CellFaceList(const t_Cell& a_Cell) {

	pCell = &a_Cell;

	if (pCell->Kind == t_CellKind::Brick) {

		for (int i = 0; i < pCell->NFaces; i++) F2V[i].setSize(4);

		// list of faces for a hexa cell
		// decomposition according cgns documentation: sids/conv.html#unst_hexa
		// vertexes (1,2,3,4,5,6,7,8) => faces 
		//(1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8)

		const t_Vert* V1 = pCell->getpVert(0);
		const t_Vert* V2 = pCell->getpVert(1);
		const t_Vert* V3 = pCell->getpVert(2);
		const t_Vert* V4 = pCell->getpVert(3);
		const t_Vert* V5 = pCell->getpVert(4);
		const t_Vert* V6 = pCell->getpVert(5);
		const t_Vert* V7 = pCell->getpVert(6);
		const t_Vert* V8 = pCell->getpVert(7);

		F2V[0][0] = V1;
		F2V[0][1] = V4;
		F2V[0][2] = V3;
		F2V[0][3] = V2;

		F2V[1][0] = V1;
		F2V[1][1] = V2;
		F2V[1][2] = V6;
		F2V[1][3] = V5;

		F2V[2][0] = V2;
		F2V[2][1] = V3;
		F2V[2][2] = V7;
		F2V[2][3] = V6;

		F2V[3][0] = V3;
		F2V[3][1] = V4;
		F2V[3][2] = V8;
		F2V[3][3] = V7;

		F2V[4][0] = V1;
		F2V[4][1] = V5;
		F2V[4][2] = V8;
		F2V[4][3] = V4;

		F2V[5][0] = V5;
		F2V[5][1] = V6;
		F2V[5][2] = V7;
		F2V[5][3] = V8;

		return;
	}

	if (pCell->Kind == t_CellKind::Tetra) {

		for (int i = 0; i < pCell->NFaces; i++) F2V[i].setSize(3);

		// list of faces for a tetra cell
		// decomposition according cgns documentation: sids/conv.html#unst_tetra
		// vertexes (1,2,3,4) => faces (1,3,2), (1,2,4), (2,3,4), (3,1,4)
		const t_Vert* V1 = pCell->getpVert(0);
		const t_Vert* V2 = pCell->getpVert(1);
		const t_Vert* V3 = pCell->getpVert(2);
		const t_Vert* V4 = pCell->getpVert(3);

		F2V[0][0] = V1;
		F2V[0][1] = V3;
		F2V[0][2] = V2;

		F2V[1][0] = V1;
		F2V[1][1] = V2;
		F2V[1][2] = V4;

		F2V[2][0] = V2;
		F2V[2][1] = V3;
		F2V[2][2] = V4;

		F2V[3][0] = V3;
		F2V[3][1] = V1;
		F2V[3][2] = V4;

		return;
	}
	hsLogMessage("t_CellFaceList: unsupported element type");
};

const t_SetOfpVerts& t_CellFaceList::getVertices(int indFace) const {

	return F2V[indFace];
};
//************************************* Face methods

t_SetOfpVerts t_Face::getVerts() const {
	t_SetOfpVerts verts;
	verts.setSize(this->NVerts);
	for (int i = 0; i < verts.size(); i++) verts[i] = pVerts[i];
	return verts;
};

void t_Zone::init_face2cell_conn(lint a_id, t_Cell& cell, int face_ind) {

	t_Face& face = Faces[a_id];

	face.Id = a_id;

	t_CellFaceList flist(cell);

	face.NVerts = flist.NVertInFace(face_ind);

	t_SetOfpVerts vertices = flist.getVertices(face_ind);

	for (int i = 0; i < face.NVerts; i++) {

		face.pVerts[i] = vertices[i];

	}


	face.pMyCell = &cell;

	face.IndLeftCellFace = face_ind;

	cell.pFaces[face_ind] = &face;

	if (cell.pCellsNeig[face_ind] != nullptr) {

		t_Cell& cell_neig = *cell.pCellsNeig[face_ind];

		face.pOppCell = &cell_neig;

		int neig_face_ind = cell.FaceIndNeig[face_ind];

		cell_neig.pFaces[neig_face_ind] = &face;

		face.IndRightCellFace = cell.FaceIndNeig[face_ind];

		face.BCKind = t_FaceBCKind::Fluid;

	}

}

//************************************* Zone methods

void t_Zone::initialize(lint a_nVerts, lint a_nCellsReal, lint a_nCellsTot) {

	nVerts = a_nVerts;
	nCellsReal = a_nCellsReal;
	nCellsTot = a_nCellsTot;

	Verts = new t_Vert[nVerts]; Cells = new t_Cell[nCellsTot];

	// zero-based ids
	for (lint i = 0; i < nVerts; i++) Verts[i].Id = i;
	for (lint i = 0; i < nCellsTot; i++) Cells[i].Id = i;
};

// make connectivity of real vertices to real vertices
// the lists of Neig Cells will be updated later with ghost cells
void t_Zone::makeVertexConnectivity() {

	t_Cell* pCell;
	t_Vert* pVert;

	// reset Neighbor counts
	for (lint i = 0; i < nVerts; i++) getpVert(i)->NNeigCells = 0;

	// first count Neighbor cells for each vertex
	for (lint i = 0; i < nCellsReal; i++) {

		pCell = getpCell(i);

		for (int j = 0; j < pCell->NVerts; j++) {

			pVert = pCell->getpVert(j);
			pVert->NNeigCells++;
		}
	}

	// allocate memory
	for (lint i = 0; i < nVerts; i++) {
		pVert = getpVert(i);
		pVert->pNeigCells = new t_Cell * [pVert->NNeigCells];
	}

	// buffer for current index in neig cell list for each vertex
	int* IndBuf = new int[nVerts];

	for (lint i = 0; i < nVerts; i++) IndBuf[i] = 0;

	// initialize reference to neighbor cells for each vertex
	for (lint i = 0; i < nCellsReal; i++) {

		pCell = getpCell(i);

		for (int j = 0; j < pCell->NVerts; j++) {

			pVert = pCell->getpVert(j);
			pVert->pNeigCells[IndBuf[pVert->Id]] = pCell;
			IndBuf[pVert->Id]++;
		}
	}

	// debug
	//hsLogMessage("Vertex 2 cell connectivity:\n");
	//for (int i = 0; i < nVerts; i++) {
	//	for (int j = 0; j < Verts[i].NNeigCells; j++)
	//		std::cout << Verts[i].pNeigCells[j]->Id << ";";
	//	std::cout << "\n";
	//}


	delete[] IndBuf;
}

// get cells that owns at least one of vertices of the particular cell face
// (excluding cell itself)
// TODO: vector push_back performance ?
// TODO: ghosts
std::vector<t_Cell*> t_Zone::getNeigCellsOfCellFace(const t_Cell& cell, int face_ind) const {

	std::vector<t_Cell*> vec_pcells(0);

	t_CellFaceList cfacelst(cell);

	t_SetOfpVerts verts = cfacelst.getVertices(face_ind);

	for (int i = 0; i < cfacelst.NVertInFace(face_ind); i++) {

		const t_Vert* Vert = verts[i];

		for (int j = 0; j < Vert->NNeigCells; j++) {

			if (Vert->pNeigCells[j]->Id != cell.Id)  vec_pcells.push_back(Vert->pNeigCells[j]);

		}

	}

	return vec_pcells;

};

int t_Zone::getFacePos(lint cell_id, const std::vector<lint> vert_ids) const {

	const t_Cell& cell = getCell(cell_id);

	t_CellFaceList flst(cell);

	t_SetOfpVerts verts; for (int i = 0; i < 4; i++) verts[i] = nullptr;
	verts.setSize(vert_ids.size());
	for (int i = 0; i < verts.size(); i++) verts[i] = getpVert(vert_ids[i]);

	for (int j = 0; j < cell.NFaces; j++) {

		t_SetOfpVerts face_verts = flst.getVertices(j);

		if (t_SetOfpVerts::cmp_weak(verts, face_verts))
			return j;

	}

	hsLogError("t_Zone:getFacePos(iCell, vert_ids) failed, icell=%ld", cell_id);
	return -1;

};

// for a given set of vertices (usually forming a face)
// find a cell that owns them all
// IMPORTANT: to be used only with abutted faces
// for inner face it will find the first cell or the second
void t_Zone::getNeigAbutCellId(const std::vector<lint>& vert_ids, lint& cell_id, int& face_pos) const {

	int npoints = vert_ids.size();

	for (int i = 0; i < npoints; i++) {

		const t_Vert& vert = getVert(vert_ids[i]);

		for (int j = 0; j < vert.NNeigCells; j++) {

			const t_Cell& cell = *vert.pNeigCells[j];

			int n_coincidence = 0;

			for (int k = 0; k < cell.NVerts; k++)
				for (int p = 0; p < npoints; p++)
					if (cell.getVert(k).Id == vert_ids[p])
						n_coincidence++;

			if (n_coincidence == npoints) {
				cell_id = cell.Id;
				face_pos = getFacePos(cell_id, vert_ids);
				return;
			}

		}

	}

	hsLogMessage("t_Zone::getNeigCellId: failed to find cell that has vertices:");
	for (int i = 0; i < npoints; i++) hsLogMessage("Vert %d:%ld", i, vert_ids[i]);

};

void t_Zone::makeCellConnectivity() {

	// 1) construct pairs of cells that have common face:
	//		a)  iterate over cells that are neighbors of face vertexes => list of "adjacent" cells
	//		b) for each of "adjacent" cells get all faces
	//		c) if face vertexes of adjacent cell coincide with vertexes of the cell, they are really adjacent

	t_Cell* pcell_base, * pcell_neig;

	t_SetOfpVerts vrtxset_base, vrtxset_neig;

	for (lint i = 0; i < nCellsReal; i++) {

		pcell_base = getpCell(i);

		t_CellFaceList cfacelst_base(*pcell_base);

		for (int j = 0; j < cfacelst_base.NFaces(); j++) {

			vrtxset_base = cfacelst_base.getVertices(j);

			std::vector<t_Cell*> vec_neig_cells = getNeigCellsOfCellFace(*pcell_base, j);

			int nneig_cells = vec_neig_cells.size();

			for (int k = 0; k < nneig_cells; k++) {

				pcell_neig = vec_neig_cells[k];

				t_CellFaceList cfacelst_neig(*pcell_neig);

				for (int p = 0; p < cfacelst_neig.NFaces(); p++) {

					t_SetOfpVerts vrtxset_neig = cfacelst_neig.getVertices(p);

					if (pcell_base->pCellsNeig[j] == nullptr) {

						// check if both faces consist of the same vertices, order not important
						if (t_SetOfpVerts::cmp_weak(vrtxset_base, vrtxset_neig)) {

							pcell_base->pCellsNeig[j] = pcell_neig;
							pcell_base->FaceIndNeig[j] = p;

							// debug messages
							//hsLogMessage("Intercell face: LeftCell_id=%d, RightCell_id=%d", pcell_base->Id, pcell_neig->Id);
							//hsLogMessage("Face Vertices: %d,%d,%d,%d", vrtxset_neig[0], vrtxset_neig[1], vrtxset_neig[2], vrtxset_neig[3]);
						}

					}

				}



			}

		}

	}

	// real-2-real connections are set, now add real-2-ghost connections
	// TODO: all zones involved here, must be a method from ghost manager...
	int ZoneID = this->idGlob;
	for (int j = 0; j < G_pMesh->nZones; j++) {

		// ghost nodes from zone j for zoneID
		const t_GhostLayer& glayer = G_GhostManager.getGhostLayer(ZoneID, j);

		if (glayer.size() > 0) {
			lint offset = G_GhostManager.calcIndOffset(ZoneID, j);
			for (int k = 0; k < glayer.size(); k++) {

				const t_Cell2GhostData gdata = glayer.data[k];

				t_Cell& cell_base = getCell(gdata.id_my);
				int face_pos_my = gdata.face_pos_my;

				t_Cell& cell_neig = *getpCell(offset + k);

				cell_base.pCellsNeig[face_pos_my] = &cell_neig;
				cell_base.FaceIndNeig[face_pos_my] = gdata.face_pos_dnr;

				// currently backward references are not required
				// TODO: this may be used later
				cell_neig.pCellsNeig[gdata.face_pos_dnr] = &cell_base;
				cell_neig.FaceIndNeig[gdata.face_pos_dnr] = face_pos_my;

			}
		}

	}

}
// make list of faces
// store only one face per intercell
// 1) count number of faces
// 2) allocate & initialize
void t_Zone::makeFaces() {

	lint nFaceMax = nCellsReal * MaxNumFacesInCell;

	bool* faces_skipped = new bool[nCellsReal * MaxNumFacesInCell];

	t_Face* faces_long = new t_Face[nCellsReal * MaxNumFacesInCell];

	for (int i = 0; i < nFaceMax; i++) faces_skipped[i] = false;

	nFaces = 0;

	// counting faces & making skip list
	for (lint i = 0; i < nCellsReal; i++) {

		const t_Cell& cell_base = getCell(i);

		t_CellFaceList cfacelst(cell_base);

		for (int j = 0; j < cell_base.NFaces; j++) {

			t_SetOfpVerts verts_base = cfacelst.getVertices(j);

			if (faces_skipped[i * MaxNumFacesInCell + j] == false) {

				nFaces++;
				// find adjacent face if this is an intercell face
				if (cell_base.pCellsNeig[j] != nullptr) {

					const t_Cell& cell_neig = *cell_base.pCellsNeig[j];
					// skip face for real cell
					if (isRealCell(cell_neig.Id)) {
						int neig_face_ind = cell_base.FaceIndNeig[j];
						faces_skipped[cell_neig.Id * MaxNumFacesInCell + neig_face_ind] = true;
					}


				}

			}

		}

	}

	// debug messages
	hsLogMessage("Zone #%d has %d faces, initializing Face list", idGlob, nFaces);

	Faces = new t_Face[nFaces];

	// initializing faces
	int iFace = 0;
	for (lint i = 0; i < nCellsReal; i++) {

		t_Cell& cell_base = getCell(i);

		for (int j = 0; j < cell_base.NFaces; j++) {

			if (faces_skipped[i * MaxNumFacesInCell + j] == false) {

				t_Face& face = Faces[iFace];

				init_face2cell_conn(iFace, cell_base, j);
				cell_base.calcFaceNormalAreaOutward(j, face.Normal, face.Area);
				iFace++;

				//debug
				if (face.BCKind == t_FaceBCKind::Fluid)
					hsLogMessage("fluid face: Face_id=%ld, owner_cell_id=%ld", iFace, cell_base.Id);

			}
		}
	}

	delete[] faces_skipped;

}

// cgns ctx prepared face_patch data
// compare all faces to faces from face patch and update their BCKind
void t_Zone::updateFacesWithBCPatch(const t_Face* face_patch, const int NFacesInPatch) {

	for (int i = 0; i < NFacesInPatch; i++) {

		const t_Face& face_in_patch = face_patch[i];

		t_SetOfpVerts vrtxset_base = face_in_patch.getVerts();

		bool face_found = false;

		//1) make list of cells that are neighbors for vertices of a current face
		for (int j = 0; j < face_in_patch.NVerts; j++) {

			const t_Vert* pVert = face_in_patch.pVerts[j];

			for (int k = 0; k < pVert->NNeigCells; k++) {
				const t_Cell* pCellNeig = pVert->pNeigCells[k];

				t_CellFaceList flist(*pCellNeig);

				for (int p = 0; p < pCellNeig->NFaces; p++) {

					t_SetOfpVerts vrtxset_neig = flist.getVertices(p);

					if (t_SetOfpVerts::cmp_weak(vrtxset_base, vrtxset_neig)) {
						// we found corresponding face
						// updating info
						pCellNeig->pFaces[p]->BCKind = face_in_patch.BCKind;
						face_found = true;
					}
				}
			}

		}
		// debug
		if (!face_found)
			hsLogMessage("Error:t_Zone::updateFacesWithBCPatch:can't find corresponding face");

	}

};

//************************************* Mesh methods

void t_Mesh::initializeFromCtx() {

	nZones = G_CGNSCtx.nZones;

	// assignZonesToProcs() will be here when multiblock is up
// now we need only G_Domain.map_iZne2cgID

	map_iZne2cgID = new int[this->nZones];

//
// Default layout: one-to-one mapping of zones indices to CGNS zone IDs
//
	for (int b = 0; b < nZones; ++b)
		map_iZne2cgID[b] = b + 1;

	Zones = new t_Zone[nZones];

	for (int i = 0; i < nZones; i++) { 

		t_Zone& zne = Zones[i];
		int cg_id = map_iZne2cgID[i];
		const t_CGNSZone& cgZne = G_CGNSCtx.cgZones[i];
		
		zne.setIdGlob(i); 

		zne.setName(cgZne.getName());


	
	}

	// bunch of code from read_cgns_mesh()
	{


		// get sizes of zones and read real cells from cgns ctx
		loadCells();

		// set up connections from verts to real cells
		makeVertexConnectivity();

		G_GhostManager.setDom(*this);
		G_GhostManager.initialize(G_CGNSCtx);

		makeCellConnectivity();

		makeFaces();


		// update mesh with bc sets
		loadBCs();

		if (checkNormalOrientations())
			hsLogMessage("check Face Normal Orientations : Ok");
		else
			hsLogMessage("Error:checkNormalOrientations failed!");

		calcUnitOstrogradResid();

		// Volume conditions info (frozen zones)
		//parseVCs(ctx);

	}



}

void t_Mesh::loadCells() {

	const t_CGNSContext& ctx = G_CGNSCtx;

	for (int iZne = 0; iZne < nZones; ++iZne)
	{
		const int& cgZneID = map_iZne2cgID[iZne];
		t_Zone& Zne = Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		cgsize_t NCellsCG = cgZne.countCells();

		// check that ncells is ok
		if (NCellsCG != cgZne.getNCells())
			hsLogError("loadCells: number of cells in cgZne is different from what was read from cgns file!");

		cgsize_t NGhosts = ctx.getNumOfGhostsForZone(cgZneID);
		// Zone stores real cells + some ghost cells

		cgsize_t NCellsTot = NCellsCG + NGhosts;

		Zne.initialize(cgZne.getNVerts(), NCellsCG, NCellsTot);

		if (NCellsCG > Zne.getnCellsTot())
			hsLogError("loadCells: failed to initialize Zne: %ld cells in CGNS Zone, %ld cells in Zone",
				NCellsCG, Zne.getnCellsTot());

		// filling in real cells
		int iCell = 0;

		for (int i = 0; i < cgZne.getSectsCell().size(); i++) {

			const t_CGSection& cg_cells = cgZne.getSectionCell(i);

			for (int j = 0; j < cg_cells.get_buf().nRows; j++) {

				t_Cell& cell = Zne.getCell(iCell);

				cell.setKind(getElementKind(cg_cells.itype));

				for (int k = 0; k < cg_cells.get_buf().nCols; k++) {
					// cgns IDs are 1-based, we use zero-based ids
					lint Vert_ID = cg_cells.get_buf().get_val(j, k) - 1;
					cell.pVerts[k] = &(Zne.getVert(Vert_ID));
				}

				iCell++;

			}
		}
		// load grid coords
		for (int i = 0; i < cgZne.getNVerts(); i++) {
			t_Vert& vert = Zne.getVert(i);
			vert.xyz[0] = cgZne.getXCoords()[i];
			vert.xyz[1] = cgZne.getYCoords()[i];
			vert.xyz[2] = cgZne.getZCoords()[i];
		}

		//std::cout << "______________________Debug, Zone Verts:\n";
		// debug output of sections
		//for (int i = 0; i < Zne.getnCells(); i++) {
		//	const t_Cell& cell = Zne.getCell(i);
		//	for (int j = 0; j < cell.NVerts; j++)
		//		std::cout << cell.getVert(j).Id<< ";";
		//	std::cout << std::endl;
		//}

	}

};

// parse BCs after they have already been read into ctx
// face lists in zones must be initialized too
void t_Mesh::loadBCs() {

	const t_CGNSContext& ctx = G_CGNSCtx;

	for (int iZne = 0; iZne < this->nZones; ++iZne)
	{
		const int cgZneID = iZne + 1;
		t_Zone& Zne = Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		for (int ipatch = 0; ipatch < cgZne.getSectsBC().size(); ipatch++) {

			const t_CGSection& fpatch_cg = cgZne.getSectionBC(ipatch);

			t_FaceBCKind bc_kind;
			bool ok = G_BCList.getBCKindBySectName(fpatch_cg.name, bc_kind);
			//ok if the patch is a bc patch (not a zone-2-zone patch)
			if (ok) {
				int nF = fpatch_cg.get_buf().nRows;
				int NVertsInFace = fpatch_cg.get_buf().nCols;
				t_Face* flist = new t_Face[nF];

				for (int iF = 0; iF < nF; iF++) {

					t_Face& face = flist[iF];

					face.BCKind = bc_kind;

					face.NVerts = NVertsInFace;

					for (int j = 0; j < NVertsInFace; j++) {

						cgsize_t iVert = fpatch_cg.get_buf().get_val(iF, j);
						// cg Id is 1-based, we make our 0-based
						face.pVerts[j] = Zne.getpVert(iVert - 1);
					}

				}

				Zne.updateFacesWithBCPatch(flist, nF);

				delete[] flist;
			}
		}
	}
};     // boundary conditions

void t_Mesh::makeVertexConnectivity() {

	for (int i = 0; i < nZones; i++) Zones[i].makeVertexConnectivity();

}


void t_Mesh::makeCellConnectivity() {

	for (int i = 0; i < nZones; i++) Zones[i].makeCellConnectivity();

}

void t_Mesh::makeFaces() {

	for (int i = 0; i < nZones; i++) Zones[i].makeFaces();

}

bool t_Mesh::checkNormalOrientations() {

	bool ok = true;

	for (int iZone = 0; iZone < this->nZones; iZone++) {

		const t_Zone& zne = this->Zones[iZone];

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			const t_Cell& cell = zne.getCell(iCell);

			t_CellFaceList flist(cell);

			for (int iFace = 0; iFace < cell.NFaces; iFace++) {

				t_Vec3 norm_cell_face;
				double area;
				cell.calcFaceNormalAreaOutward(iFace, norm_cell_face, area);

				// compare computed normal & area to face normal & area
				const t_Face& face = *cell.pFaces[iFace];
				double scal_prod = norm_cell_face.dot(face.Normal);
				if (&cell == face.pMyCell) {
					ok = ok && (scal_prod > 0.999) && (scal_prod < 1.001);
					//hsLogMessage("Face norm & cell face norm must be the same (scal prod=+1), computed:%lf", scal_prod);
				}
				else if (&cell == face.pOppCell) {
					ok = ok && (scal_prod > -1.001) && (scal_prod < -0.999);
					//hsLogMessage("Face norm & cell face norm must be opposite (scal prod=-1), computed:%lf", scal_prod);
				}
				else { hsLogMessage("Error:checkNormalOrientations():Broken Face2Cell connectivity"); }
			}
		}

	}

	return ok;

}
// calculate residual for a Ostrogradsky theorem
// for all cells
double t_Mesh::calcUnitOstrogradResid() {

	double max_resid = 0.0;
	double resid;
	double area;
	t_Vec3 norm, sum;
	for (int iZone = 0; iZone < this->nZones; iZone++) {

		const t_Zone& zne = this->Zones[iZone];

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			const t_Cell& cell = zne.getCell(iCell);

			sum.set(0, 0, 0);

			for (int k = 0; k < cell.NFaces; k++) {

				norm = cell.getFaceNormalOutward(k);
				area = cell.getFace(k).Area;
				sum += norm * area;
			}

			resid = sum.norm();
			if (resid > max_resid) max_resid = resid;
		}

	}
	hsLogMessage("Check volumes (Ostrogradsky): Max resid = %lf", max_resid);
	return max_resid;
}



