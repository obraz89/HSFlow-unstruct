#include "logging.h"
#include "common_data.h"
#include "common_procs.h"

// mainly for debugging
#include <iostream>

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

t_CellEdgeList::t_CellEdgeList(const t_Cell& a_Cell) :Cell(a_Cell) {

	if (Cell.Kind == t_CellKind::Brick) {

		// list of edges for a hexa cell
		// decomposition according cgns documentation: sids/conv.html#unst_hexa
		// vertexes (1,2,3,4,5,6,7,8) => edges 
		// (1,2), (2,3), (3,4), (4,1),
		// (1,5), (2,6), (3,7), (4,8),
		// (5,6), (6,7), (7,8), (8,5)

		const lint& V1 = Cell.getVert(0).Id;
		const lint& V2 = Cell.getVert(1).Id;
		const lint& V3 = Cell.getVert(2).Id;
		const lint& V4 = Cell.getVert(3).Id;
		const lint& V5 = Cell.getVert(4).Id;
		const lint& V6 = Cell.getVert(5).Id;
		const lint& V7 = Cell.getVert(6).Id;
		const lint& V8 = Cell.getVert(7).Id;

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

	if (Cell.Kind == t_CellKind::Tetra) {

		// list of edges for a tetra cell
		// decomposition according cgns documentation: sids/conv.html#unst_tetra
		// vertexes (1,2,3,4) => edges (1,2),(2,3),(3,1),(1,4),(2,4),(3,4)

		const lint& V1 = Cell.getVert(0).Id;
		const lint& V2 = Cell.getVert(1).Id;
		const lint& V3 = Cell.getVert(2).Id;
		const lint& V4 = Cell.getVert(3).Id;

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

t_CellFaceList::t_CellFaceList(const t_Cell& a_Cell) :Cell(a_Cell) {

	if (Cell.Kind == t_CellKind::Brick) {

		for (int i = 0; i < Cell.NFaces; i++) ArrNumOfVertsInFaces[i] = 4;

		// list of faces for a hexa cell
		// decomposition according cgns documentation: sids/conv.html#unst_hexa
		// vertexes (1,2,3,4,5,6,7,8) => faces 
		//(1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8)

		const lint& V1 = Cell.getVert(0).Id;
		const lint& V2 = Cell.getVert(1).Id;
		const lint& V3 = Cell.getVert(2).Id;
		const lint& V4 = Cell.getVert(3).Id;
		const lint& V5 = Cell.getVert(4).Id;
		const lint& V6 = Cell.getVert(5).Id;
		const lint& V7 = Cell.getVert(6).Id;
		const lint& V8 = Cell.getVert(7).Id;

		F2V.get_val(0, 0) = V1;
		F2V.get_val(0, 1) = V4;
		F2V.get_val(0, 2) = V3;
		F2V.get_val(0, 3) = V2;

		F2V.get_val(1, 0) = V1;
		F2V.get_val(1, 1) = V2;
		F2V.get_val(1, 2) = V6;
		F2V.get_val(1, 3) = V5;

		F2V.get_val(2, 0) = V2;
		F2V.get_val(2, 1) = V3;
		F2V.get_val(2, 2) = V7;
		F2V.get_val(2, 3) = V6;

		F2V.get_val(3, 0) = V3;
		F2V.get_val(3, 1) = V4;
		F2V.get_val(3, 2) = V8;
		F2V.get_val(3, 3) = V7;

		F2V.get_val(4, 0) = V1;
		F2V.get_val(4, 1) = V5;
		F2V.get_val(4, 2) = V8;
		F2V.get_val(4, 3) = V4;

		F2V.get_val(5, 0) = V5;
		F2V.get_val(5, 1) = V6;
		F2V.get_val(5, 2) = V7;
		F2V.get_val(5, 3) = V8;

		return;
	}

	if (Cell.Kind == t_CellKind::Tetra) {

		for (int i = 0; i < Cell.NFaces; i++) ArrNumOfVertsInFaces[i] = 3;

		// list of faces for a tetra cell
		// decomposition according cgns documentation: sids/conv.html#unst_tetra
		// vertexes (1,2,3,4) => faces (1,2,3), (1,2,4), (2,3,4), (3,1,4)
		const lint& V1 = Cell.getVert(0).Id;
		const lint& V2 = Cell.getVert(1).Id;
		const lint& V3 = Cell.getVert(2).Id;
		const lint& V4 = Cell.getVert(0).Id;

		F2V.get_val(0, 0) = V1;
		F2V.get_val(0, 1) = V2;
		F2V.get_val(0, 2) = V3;

		F2V.get_val(1, 0) = V1;
		F2V.get_val(1, 1) = V2;
		F2V.get_val(1, 2) = V4;

		F2V.get_val(2, 0) = V2;
		F2V.get_val(2, 1) = V3;
		F2V.get_val(2, 2) = V4;

		F2V.get_val(3, 0) = V3;
		F2V.get_val(3, 1) = V1;
		F2V.get_val(3, 2) = V4;

		return;
	}
	hsLogMessage("t_CellFaceList: unsupported element type");
};

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
	for (int i = 0; i < nVerts; i++) getpVert(i)->NNeigCells = 0;

	// first count Neighbor cells for each vertex
	for (int i = 0; i < nCells; i++) {

		pCell = getpCell(i);

		for (int j = 0; j < pCell->NVerts; j++) {

			pVert = pCell->getpVert(j);
			pVert->NNeigCells++;
		}
	}

	// allocate memory
	for (int i = 0; i < nVerts; i++) {
		pVert = getpVert(i);
		pVert->pNeigCells = new t_Cell*[pVert->NNeigCells];
	}

	// buffer for current index in neig cell list for each vertex
	int* IndBuf = new int[nVerts];

	for (int i = 0; i < nVerts; i++) IndBuf[i] = 0;

	// initialize reference to neighbor cells for each vertex
	for (int i = 0; i < nCells; i++) {

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

void t_Zone::getNeigCellsOfVertices() {



};

void t_Zone::makeFaceList() {

	// make the list of all faces in each zone
	// also make lists for bcs
	// bcs are detected as faces that have only 1 cell attached to it

	// the major steps are:

	// 1) construct list of neighbors for each vertex

	getNeigCellsOfVertices();

	// 2) construct pairs of cells that have common face:
	//		a)  iterate over vertexes that are neighbors of face vertexes => list of "adjacent" cells
	//		b) for each of "adjacent" cells get all faces
	//		c) if face vertexes of adjacent cell coincide with vertexes of the cell, they really are adjacent


}

void t_Domain::makeVertexConnectivity() {

	for (int i = 0; i < nZones; i++) Zones[i].makeVertexConnectivity();

}


void t_Domain::makeFaceLists() {

	for (int i = 0; i < nZones; i++) Zones[i].makeFaceList();

}

