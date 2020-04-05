#pragma once

#include <map>

#include "cgnslib.h"

#include "common_data.h"
#include "common_procs.h"

using t_BufCGSize = t_BufInds<cgsize_t>;

// array of typical cgns elements like face, cell
// or smth else defined as set of vertices
struct t_CGElemArray {

	CG_ElementType_t itype;

	// container to store indices of vertices
	t_BufCGSize _buf;
	// array to store cg ids of elements
	cgsize_t* _ids;

	t_BufCGSize& get_buf() { return _buf; }
	const t_BufCGSize& get_buf() const{ return _buf; }

	cgsize_t* get_buf_data() { return _buf.data(); }
	const cgsize_t* get_buf_data() const{ return _buf.data(); }

	void alloc(cgsize_t n_elems, int n_verts_in_elem) {

		_buf.allocate(n_elems, n_verts_in_elem);
		delete[] _ids;
		_ids = new cgsize_t[n_elems];
	}

	t_CGElemArray() :_buf(0, 0), itype(CG_ElementTypeNull) { _ids = new cgsize_t[0]; }

	~t_CGElemArray() { delete[] _ids; }


};

// patch of fixed element type, mixed elements not supported
struct t_CGFacePatch : public t_CGElemArray
{

	std::string nameOfSect;

	t_CGFacePatch(const std::string& name):nameOfSect(name), t_CGElemArray() { ; }
};

struct t_CGNSZone
{
	// vector of arrays of cells {[Tetra], [HEXA-8], ...}
	std::vector<t_CGElemArray*> pCellSets;

	// bc patches
	std::vector<t_CGFacePatch*> pPatchesBC;
	// abutting patches
	std::vector<t_CGFacePatch*> pPatchesAbut;

	cgsize_t countCells() {
		cgsize_t N=0;
		for (int i = 0; i < pCellSets.size(); i++) 
			N += pCellSets[i]->get_buf().nRows;
		return N;
	}

	t_CGNSZone() : pPatchesAbut(), pPatchesBC(), pCellSets() { ; }
	~t_CGNSZone() { 
		for (int i = 0; i < pCellSets.size(); i++) delete pCellSets[i];
		for (int i = 0; i < pPatchesAbut.size(); i++) delete pPatchesAbut[i]; 
		for (int i = 0; i < pPatchesBC.size(); i++) delete pPatchesBC[i];
	}
};

struct t_CGNSContext
{
	int iFile, iBase;
	std::map<std::string, int>  map_ZneName2Idx;

	t_CGNSZone* cgZones;

	t_CGNSContext() : cgZones(NULL) { ; }
	~t_CGNSContext() { delete[] cgZones; }
};

int read_cgns_mesh();