#include "logging.h"
#include "common_data.h"
#include "common_procs.h"

// mainly for debugging
#include <iostream>

#include <vector>

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

void t_CellEdgeList::init(const t_Cell& a_Cell){

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

t_CellFaceList::t_CellFaceList(const t_Cell& a_Cell){

	pCell = &a_Cell;

	if (pCell->Kind == t_CellKind::Brick) {

		for (int i = 0; i < pCell->NFaces; i++) F2V[i].setSize(4);

		// list of faces for a hexa cell
		// decomposition according cgns documentation: sids/conv.html#unst_hexa
		// vertexes (1,2,3,4,5,6,7,8) => faces 
		//(1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8)

		const lint& V1 = pCell->getVert(0).Id;
		const lint& V2 = pCell->getVert(1).Id;
		const lint& V3 = pCell->getVert(2).Id;
		const lint& V4 = pCell->getVert(3).Id;
		const lint& V5 = pCell->getVert(4).Id;
		const lint& V6 = pCell->getVert(5).Id;
		const lint& V7 = pCell->getVert(6).Id;
		const lint& V8 = pCell->getVert(7).Id;

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
		// vertexes (1,2,3,4) => faces (1,2,3), (1,2,4), (2,3,4), (3,1,4)
		const lint& V1 = pCell->getVert(0).Id;
		const lint& V2 = pCell->getVert(1).Id;
		const lint& V3 = pCell->getVert(2).Id;
		const lint& V4 = pCell->getVert(3).Id;

		F2V[0][0] = V1;
		F2V[0][1] = V2;
		F2V[0][2] = V3;

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

const t_SetIndF2V& t_CellFaceList::getVertices(int indFace) const{

	return F2V[indFace];
};
//************************************* Face methods

void t_Zone::init_face2cell_conn(lint a_id, t_Cell& cell, int face_ind) {
	
	t_Face& face = Faces[a_id];

	face.Id = a_id;

	t_CellFaceList flist(cell);

	face.NVerts = flist.NVertInFace(face_ind);

	t_SetIndF2V vertices = flist.getVertices(face_ind);

	for (int i = 0; i < face.NVerts; i++) {

		face.pVerts[i] = getpVert(vertices[i]);

	}

	// TODO: for now i leave left cell as base cell
	// and right cell as neig cell
	// so left cell is always defined but right cell
	// can be missing (BC face)
	face.pLeftCell = &cell;

	face.IndLeftCellFace = face_ind;

	cell.pFaces[face_ind] = &face;

	if (cell.pCellsNeig[face_ind] != nullptr) {

		t_Cell& cell_neig = *cell.pCellsNeig[face_ind];

		face.pRightCell = &cell_neig;

		int neig_face_ind = cell.FaceIndNeig[face_ind];

		cell_neig.pFaces[neig_face_ind] = &face;

		face.IndRightCellFace = cell.FaceIndNeig[face_ind];

		face.BCKind = t_FaceBCKind::Fluid;

	}

}

//************************************* Zone methods

void t_Zone::initialize(lint a_nVerts, lint a_nCells) {

	nVerts = a_nVerts;
	nCells = a_nCells;

	Verts = new t_Vert[nVerts]; Cells = new t_Cell[nCells];

	// zero-based ids
	for (lint i = 0; i < nVerts; i++) Verts[i].Id = i;
	for (lint i = 0; i < nCells; i++) Cells[i].Id = i;
};

void t_Zone::makeVertexConnectivity() {

	t_Cell* pCell;
	t_Vert* pVert;

	// reset Neighbor counts
	for (lint i = 0; i < nVerts; i++) getpVert(i)->NNeigCells = 0;

	// first count Neighbor cells for each vertex
	for (lint i = 0; i < nCells; i++) {

		pCell = getpCell(i);

		for (int j = 0; j < pCell->NVerts; j++) {

			pVert = pCell->getpVert(j);
			pVert->NNeigCells++;
		}
	}

	// allocate memory
	for (lint i = 0; i < nVerts; i++) {
		pVert = getpVert(i);
		pVert->pNeigCells = new t_Cell*[pVert->NNeigCells];
	}

	// buffer for current index in neig cell list for each vertex
	int* IndBuf = new int[nVerts];

	for (lint i = 0; i < nVerts; i++) IndBuf[i] = 0;

	// initialize reference to neighbor cells for each vertex
	for (lint i = 0; i < nCells; i++) {

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
std::vector<t_Cell*> t_Zone::getNeigCellsOfCellFace(const t_Cell& cell, int face_ind) const{

	std::vector<t_Cell*> vec_pcells(0);

	t_CellFaceList cfacelst(cell);

	t_SetIndF2V verts_ids = cfacelst.getVertices(face_ind);

	for (int i = 0; i < cfacelst.NVertInFace(face_ind); i++) {

		const t_Vert& Vert = getVert(verts_ids[i]);

		for (int j = 0; j < Vert.NNeigCells; j++) {

			if (Vert.pNeigCells[j]->Id != cell.Id)  vec_pcells.push_back(Vert.pNeigCells[j]);

		}

	}

	return vec_pcells;

};

void t_Zone::makeCellConnectivity() {

	// make the list of all faces in each zone
	// also detect bc faces : faces that have only 1 cell attached to it

	// 1) construct pairs of cells that have common face:
	//		a)  iterate over cells that are neighbors of face vertexes => list of "adjacent" cells
	//		b) for each of "adjacent" cells get all faces
	//		c) if face vertexes of adjacent cell coincide with vertexes of the cell, they are really adjacent

	// count number of faces
	lint nFaces = 0;

	t_Cell *pcell_base, *pcell_neig;

	t_SetIndF2V vrtxset_base, vrtxset_neig;

	for (lint i = 0; i < nCells; i++) {

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

					t_SetIndF2V vrtxset_neig = cfacelst_neig.getVertices(p);

					if (pcell_base->pCellsNeig[j] == nullptr) {

						// check if both faces consist of the same vertices, order not important
						if (t_SetIndF2V::cmp_weak(vrtxset_base, vrtxset_neig)) {

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

}
// make list of faces
// store only one face per intercell
// 1) count number of faces
// 2) allocate & initialize
void t_Zone::makeFaces() {

	lint nFaceMax = nCells*MaxNumFacesInCell;

	bool* faces_skipped = new bool[nCells*MaxNumFacesInCell];

	t_Face* faces_long = new t_Face[nCells*MaxNumFacesInCell];

	for (int i = 0; i < nFaceMax; i++) faces_skipped[i] = false;

	nFaces = 0;

	// counting faces & making skip list
	for (lint i = 0; i < nCells; i++) {

		const t_Cell& cell_base = getCell(i);

		t_CellFaceList cfacelst(cell_base);

		for (int j = 0; j < cell_base.NFaces; j++) {

			t_SetIndF2V verts_base = cfacelst.getVertices(j);

			if (faces_skipped[i*MaxNumFacesInCell + j] == false) {

				nFaces++;
				// find adjacent face if this is an intercell face
				if (cell_base.pCellsNeig[j] != nullptr) {

					const t_Cell& cell_neig = *cell_base.pCellsNeig[j];
					int neig_face_ind = cell_base.FaceIndNeig[j];
					faces_skipped[cell_neig.Id*MaxNumFacesInCell + neig_face_ind] = true;

				}

			}

		}

	}

	// debug messages
	//hsLogMessage("Zone has %d faces", nFaces);

	Faces = new t_Face[nFaces];

	// initializing faces
	int iFace = 0;
	for (lint i = 0; i < nCells; i++) {

		t_Cell& cell_base = getCell(i);

		for (int j = 0; j < cell_base.NFaces; j++) {

			if (faces_skipped[i*MaxNumFacesInCell + j] == false) {

				t_Face& face = Faces[iFace];

				init_face2cell_conn(iFace, cell_base, j);
				iFace++;

			}
		}
	}

	delete[] faces_skipped;

}

void t_Domain::makeVertexConnectivity() {

	for (int i = 0; i < nZones; i++) Zones[i].makeVertexConnectivity();

}


void t_Domain::makeCellConnectivity() {

	for (int i = 0; i < nZones; i++) Zones[i].makeCellConnectivity();

}

void t_Domain::makeFaces() {

	for (int i = 0; i < nZones; i++) Zones[i].makeFaces();

}

