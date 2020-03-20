#include "logging.h"
#include "common_data.h"
#include "common_procs.h"

void t_Zone::initialize(lint a_nVerts, lint a_nCells) {

	nVerts = a_nVerts;
	nCells = a_nCells;

	Verts = new t_Vert[nVerts]; Cells = new t_Cell[nCells];

	// for now keep ids from cgns, they are 1-based
	for (lint i = 1; i <= nVerts; i++) Verts[i].Id = i;
	for (lint i = 1; i <= nCells; i++) Cells[i].Id = i;
};



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

void t_Domain::makeFaceLists() {

	for (int i = 0; i < nZones; i++) Zones[i].makeFaceList();

}

// list of faces for a tetra cell
// decomposition according cgns documentation: sids/conv.html#unst_tetra
// vertexes (1,2,3,4) => faces (1,2,3), (1,2,4), (2,3,4), (3,1,4)
void getFacesOfTetra(const lint* const Vertexes, t_BufIndsStat<lint, 4, 3>& faces) {

	const lint& V1 = Vertexes[0];
	const lint& V2 = Vertexes[1];
	const lint& V3 = Vertexes[2];
	const lint& V4 = Vertexes[3];

	faces.get_val(0, 0) = V1;
	faces.get_val(0, 1) = V2;
	faces.get_val(0, 2) = V3;

	faces.get_val(1, 0) = V1;
	faces.get_val(1, 1) = V2;
	faces.get_val(1, 2) = V4;

	faces.get_val(2, 0) = V2;
	faces.get_val(2, 1) = V3;
	faces.get_val(2, 2) = V4;

	faces.get_val(3, 0) = V3;
	faces.get_val(3, 1) = V1;
	faces.get_val(3, 2) = V4;

}

// list of faces for a hexa cell
// decomposition according cgns documentation: sids/conv.html#unst_hexa
// vertexes (1,2,3,4,5,6,7,8) => faces 
//(1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8)
void getFacesOfHexa(const lint* const Vertexes, t_BufIndsStat<lint, 6, 4>& faces) {

	const lint& V1 = Vertexes[0];
	const lint& V2 = Vertexes[1];
	const lint& V3 = Vertexes[2];
	const lint& V4 = Vertexes[3];
	const lint& V5 = Vertexes[4];
	const lint& V6 = Vertexes[5];
	const lint& V7 = Vertexes[6];
	const lint& V8 = Vertexes[7];

	faces.get_val(0, 0) = V1;
	faces.get_val(0, 1) = V4;
	faces.get_val(0, 2) = V3;
	faces.get_val(0, 3) = V2;

	faces.get_val(1, 0) = V1;
	faces.get_val(1, 1) = V2;
	faces.get_val(1, 2) = V6;
	faces.get_val(1, 3) = V5;

	faces.get_val(2, 0) = V2;
	faces.get_val(2, 1) = V3;
	faces.get_val(2, 2) = V7;
	faces.get_val(2, 3) = V6;

	faces.get_val(3, 0) = V3;
	faces.get_val(3, 1) = V4;
	faces.get_val(3, 2) = V8;
	faces.get_val(3, 3) = V7;

	faces.get_val(4, 0) = V1;
	faces.get_val(4, 1) = V5;
	faces.get_val(4, 2) = V8;
	faces.get_val(4, 3) = V4;

	faces.get_val(5, 0) = V5;
	faces.get_val(5, 1) = V6;
	faces.get_val(5, 2) = V7;
	faces.get_val(5, 3) = V8;

}