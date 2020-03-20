#pragma once

#include <map>

#include "cgnslib.h"

class t_BufCGSize {
private:
	// TODO: for some reason, vs doesnt remove default constructor with
	// t_BufCGSize() = delete;
	// making it inaccessible old style
	t_BufCGSize() {};

	cgsize_t* buf;
public:
	cgsize_t nRows, nCols;
	cgsize_t* data() { return buf; }
	//t_BufCGSize() = delete;
	t_BufCGSize(t_BufCGSize&) = delete;
	t_BufCGSize(cgsize_t a_NR, cgsize_t a_NC) { nRows = a_NR; nCols = a_NC; buf = new cgsize_t[nRows*nCols]; }
	void allocate(cgsize_t a_NR, cgsize_t a_NC) { delete[] buf;  nRows = a_NR; nCols = a_NC; buf = new cgsize_t[nRows*nCols];};
	cgsize_t get_val(cgsize_t i, cgsize_t j) { return *(buf + i*nCols + j); };
	~t_BufCGSize() { delete[] buf; }

};

struct t_CGNSZone
{

	t_BufCGSize cells;


	// Abutted faces data. The face consists of one or many patches
	struct t_FacePatch
	{
		//

		// Ctor
		t_FacePatch() { ; }
	};

	// Zone connectivity info
	int nPatches;
	t_FacePatch* Patches;

	t_CGNSZone() : nPatches(0), Patches(NULL), cells(0,0) { ; }
	~t_CGNSZone() { delete[] Patches; }
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