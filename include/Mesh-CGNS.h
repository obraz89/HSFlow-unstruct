#pragma once

#include <map>

struct t_CGNSZone
{

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

	t_CGNSZone() : nPatches(0), Patches(NULL) { ; }
	~t_CGNSZone() { delete[] Patches; }
};

struct t_CGNSContext
{
	int fileID, iBase;
	std::map<std::string, int>  map_zoneName_id;

	t_CGNSZone* cgZones;

	t_CGNSContext() : cgZones(NULL) { ; }
	~t_CGNSContext() { delete[] cgZones; }
};

int read_cgns_mesh();