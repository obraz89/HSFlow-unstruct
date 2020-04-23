/**************************************
// Name:        ghost-common.h
// Purpose:     interface of basic ghost manager
// Author:      A. Obraz
**************************************/
/**************************************
Ghost manager handles connections
between zones by layers of ghost cells
**************************************/
#pragma once

#include "common_data.h"

#include "mesh.h"

#include "CGNS-ctx.h"

// sell-2-ghost connection data
// real cell C from some zone (Z_my) abutting ghost cell G via face F
// id_my - id of cell C in this zone (Z_my)  
// id_dnr - id of cell G in donor zone (Z_dnr)
// face_pos_my - index of face position for face F in cell C
// NB: ids of zones (Z_my, Z_dnr) are stored automatically in Ghost Manager
struct t_Cell2GhostData {

	lint id_my;
	lint id_dnr;
	int face_pos_my;
	int face_pos_dnr;

};

struct t_GhostLayer {

	std::vector<t_Cell2GhostData> data;
	
	cgsize_t size() const{ 
		return data.size(); 
	}

	void resize(lint i) {
		data.resize(i);
	}

};

// Zones are connected via layers of ghost cells;
// t_GhostManager stores connection of ghosts to real nodes

// this is a single-proc variant
// TODO: this can be upgraded to MPI version

// In MPI version, data should be serialized
// in single-proc version, we can pass complex data types
class t_GhostMngBase {

protected:

	t_Mesh* _pDom;

	// connections between zones
	// in MPI case, should be initialized via receive from master
	// stores everyone-2-everyone 
	// not a very big overhead as most of connections are empty
	std::vector<t_GhostLayer*> _pGLayers;


public:

	int getPlainInd(int i, int j) const{ return i * _pDom->nZones + j; }

	const t_GhostLayer& getGhostLayer(int i, int j) const{ 

#ifdef _DEBUG
		if ( i<0 || j<0 || i>_pDom->nZones-1 || j>_pDom->nZones - 1)
			hsLogError("t_GhostManager:getGhostLayer(i,j): bad indices: i=%d, j=%d",i,j);
#endif

		return *_pGLayers[getPlainInd(i, j)]; 
	}

	lint calcNumOfGhostCells(int zone_id) const;

	lint calcIndOffset(int zone_id_my, int zone_id_dnr) const;

	void initialize(const t_CGNSContext& ctx);

	void getGhostsZiFromZj_Neig(const t_CGNSContext& ctx, int cgZneID_I, int cgZneID_J,
		t_GhostLayer& glayer) const;

	virtual void setDom(t_Mesh& a_dom) {

		_pDom = &a_dom;

		_pGLayers.resize(_pDom->nZones * _pDom->nZones, nullptr);

	}

	t_GhostMngBase(): _pGLayers(), _pDom(nullptr) {}

	// 
	virtual void exchangeGeomData();

	// physical methods implemented by children
	virtual void exchangeCSV() = 0;

	virtual ~t_GhostMngBase() { for (int i = 0; i < _pGLayers.size(); i++) delete _pGLayers[i]; }
};

extern t_GhostMngBase* G_pGhostMngBase ;