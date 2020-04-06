#include "ghost_manager.h"

t_GhostManager G_GhostManager;

void t_GhostManager::initialize(const t_CGNSContext& ctx) {

	for (int i = 0; i < _pDom->nZones; i++) {

		for (int j = 0; j < _pDom->nZones; j++) {

			int cgZneID_I = i + 1;
			int cgZneID_J = j + 1;

			t_GhostLayer* glayer = new t_GhostLayer();

			ctx.getGhostsZiFromZj_Neig(cgZneID_I, cgZneID_J, glayer->ids_my, glayer->ids_dnr);

			_pGLayers[getPlainInd(i, j)] = glayer;

		}

	}

};

lint t_GhostManager::calcNumOfGhostCells(int zone_id) const{

	lint num_of_ghosts = 0;

	for (int j = 0; j < -_pDom->nZones; j++) {

		t_GhostLayer* glayer = _pGLayers[getPlainInd(zone_id, j)];
		num_of_ghosts += glayer->ids_dnr.size();

	}

	return num_of_ghosts;

};