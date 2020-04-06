#pragma once

#include "common_data.h"

#include "Mesh-CGNS.h"

// ghost cells for some zone (My) from zone (Dnr)
// ids_dnr - list of donor cells ids
// ids_my - list of abutting cells from this zone
// 
struct t_GhostLayer {

	std::vector<lint> ids_my;
	std::vector<lint> ids_dnr;

};

// Zones are connected via layers of ghost cells;
// t_GhostManager stores connection of ghosts to real nodes
// and also does exchanges of information between zones

// this is a single-proc variant
// TODO: this can be upgraded to MPI version

// In MPI version, data should be serialized
// in single-proc version, we can pass complex data types
class t_GhostManager {

	t_Domain* _pDom;

	// list of zones to manage
	std::vector<t_Zone*> _pMyZones;

	// connections between zones
	// in MPI case, should be initialized via receive from master
	// stores everyone-2-everyone 
	// not a very big overhead as most of connections are empty
	std::vector<t_GhostLayer*> _pGLayers;


public:

	int getPlainInd(int i, int j) const{ return i * _pDom->nZones + j; }

	lint calcNumOfGhostCells(int zone_id) const;

	void initialize(const t_CGNSContext& ctx);

	void setDom(t_Domain& a_dom) {

		_pDom = &a_dom;

		_pGLayers.resize(_pDom->nZones * _pDom->nZones, nullptr);

		// single proc
		for (int i = 0; i < _pDom->nZones; i++)
			_pMyZones.push_back(&_pDom->Zones[i]);
	}


	t_GhostManager():_pMyZones(0), _pGLayers() {}



	~t_GhostManager() { for (int i = 0; i < _pGLayers.size(); i++) delete _pGLayers[i]; }
};

extern t_GhostManager G_GhostManager ;