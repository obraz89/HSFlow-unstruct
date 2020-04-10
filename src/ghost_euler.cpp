#include "ghost_euler.h"

t_GhostMngEuler G_GhostMngEu;

// this is like send and receive 
void  t_GhostMngEuler::exchangeCSV() {

	for (int i = 0; i < _pDom->nZones; i++) {

		for (int j = 0; j < _pDom->nZones; j++) {

			t_Zone& ZoneMy = _pDom->Zones[i];
			t_Zone& ZoneDnr = _pDom->Zones[j];

			t_GhostLayer* glayer = _pGLayers[getPlainInd(i, j)];

			for (int k = 0; k < glayer->size(); k++) {

				lint offset = calcIndOffset(i, j);
				_pDomEu->getCellCSV(i, offset + k) = _pDomEu->getCellCSV(j, glayer->data[k].id_dnr);
				//ZoneMy.getCell(offset + k).ConsVars = ZoneDnr.getCell(glayer->data[k].id_dnr).ConsVars;

			}

		}
	}

};