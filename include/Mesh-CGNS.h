#pragma once

#include <map>

#include "cgnslib.h"

#include "common_data.h"
#include "common_procs.h"

using t_BufCGSize = t_BufInds<cgsize_t>;

struct t_CGNSZone
{

	// TODO: later there will be mixed elements, and sections must be joint together
	// for now assuming one type of elements read from one section
	t_BufCGSize cells;

	CGNS_ENUMT(ElementType_t) itype;


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