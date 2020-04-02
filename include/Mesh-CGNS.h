#pragma once

#include <map>

#include "cgnslib.h"

#include "common_data.h"
#include "common_procs.h"

using t_BufCGSize = t_BufInds<cgsize_t>;

// patch of fixed element type, mixed elements not supported
struct t_CGFacePatch
{

	std::string nameOfSect;

	CG_ElementType_t itype;

	t_BufCGSize _data;

	cgsize_t* data() { return _data.data(); }

	t_CGFacePatch(const std::string& name):nameOfSect(name), _data(0,0) { ; }
	// TODO: compiler wants it, don't know why
	t_CGFacePatch() :_data(0, 0), nameOfSect("") { }
	t_CGFacePatch(t_CGFacePatch& p):_data(0,0){ nameOfSect = p.nameOfSect; }
};

struct t_CGNSZone
{

	// TODO: later there will be mixed elements, and sections must be joint together
	// for now assuming one type of elements read from one section
	t_BufCGSize cells;

	CGNS_ENUMT(ElementType_t) itype;

	// Zone connectivity info
	std::vector<t_CGFacePatch*> pFacePatches;

	t_CGNSZone() : pFacePatches(), cells(0,0) { ; }
	~t_CGNSZone() { for (int i = 0; i < pFacePatches.size(); i++) delete pFacePatches[i]; }
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