/**************************************
// Name:        ghost-common.cpp
// Purpose:     Basics of ghost manager
// Author:      A. Obraz
**************************************/
/**************************************
Ghost manager handles connections 
between zones by layers of ghost cells
**************************************/
#include "ghost_common.h"

#include "CGNS-ctx.h"

t_GhostMngBase* G_pGhostMngBase;

// initialize directly from cgns ctx
// at this moment mesh is not fully initialized
void t_GhostMngBase::initialize(const t_CGNSContext& ctx) {

	if (_pDom->nZones == 0)
		hsLogMessage("Error: Initializing ghost manager while domain is not initialized!");

	_pGLayers.resize(_pDom->nZones * _pDom->nZones, nullptr);

	for (int i = 0; i < _pDom->nZones; i++) {

		for (int j = 0; j < _pDom->nZones; j++) {

			int cgZneID_I = i + 1;
			int cgZneID_J = j + 1;

			t_GhostLayer* glayer = new t_GhostLayer();

			getGhostsZiFromZj_Neig(ctx, cgZneID_I, cgZneID_J, *glayer);

			_pGLayers[getPlainInd(i, j)] = glayer;

		}

	}

};

// converting elem2elem connectivity (which is face2face)
// into cell2cell connectivity;
// ids_dnr : store ids of cells for zone_J that are neighbors of abutting cells from zone_I (ghosts)
// ids_my : store ids of cells for zone_I that are abutting with ghosts
// as a result zone_I.getCell(ids_my[k]) and zone_J.getCell(ids_dnr[k]) are abutting
// first one is real for zone_I, second is ghost for zone_I
// IMPORTANT: this function make calls to t_Zone functions
// that require Vertex connectivity being initialized
void t_GhostMngBase::getGhostsZiFromZj_Neig(const t_CGNSContext& ctx, int cgZneID_I, int cgZneID_J,
	t_GhostLayer& glayer) const {

	const t_CGNSZone& cgZneMy = ctx.cgZones[cgZneID_I - 1];
	const t_CGNSZone& cgZneDnr = ctx.cgZones[cgZneID_J - 1];

	const std::vector<t_CGConnSet*>& pConns = cgZneMy.getConns();

	glayer.resize(0);

	for (int i = 0; i < pConns.size(); i++) {

		const t_CGConnSet& conn = *pConns[i];

		if (conn.getZoneIDDnr() == cgZneID_J) {

			for (int j = 0; j < conn.get_buf().nCols; j++) {
				cgsize_t id_my, id_dnr;
				conn.getConnIds(j, id_my, id_dnr);

				// now we need to find cell that has that face...

				// first get ids of nodes
				std::vector<cgsize_t> verts_cg_ids_my = cgZneMy.getVertsOfElem(id_my);
				std::vector<cgsize_t> verts_cg_ids_dnr = cgZneDnr.getVertsOfElem(id_dnr);

				// get zero-based indices of vertices
				std::vector<cgsize_t> verts_ids_my = verts_cg_ids_my;
				for (int k = 0; k < verts_ids_my.size(); k++) verts_ids_my[k] -= 1;

				std::vector<cgsize_t> verts_ids_dnr = verts_cg_ids_dnr;
				for (int k = 0; k < verts_ids_dnr.size(); k++) verts_ids_dnr[k] -= 1;

				// now we need to use vertexConnectivity
				const t_Zone& ZneMy = G_pDom->Zones[cgZneID_I - 1];
				const t_Zone& ZneDnr = G_pDom->Zones[cgZneID_J - 1];

				int face_pos_my, face_pos_dnr;

				lint cell_id_my, cell_id_dnr;

				ZneMy.getNeigAbutCellId(verts_ids_my, cell_id_my, face_pos_my);
				ZneDnr.getNeigAbutCellId(verts_ids_dnr, cell_id_dnr, face_pos_dnr);

				t_Cell2GhostData data_point;

				data_point.id_my = cell_id_my;
				data_point.id_dnr = cell_id_dnr;
				data_point.face_pos_my = face_pos_my;
				data_point.face_pos_dnr = face_pos_dnr;


				glayer.data.push_back(data_point);

			}


		}

	}

};

lint t_GhostMngBase::calcNumOfGhostCells(int zone_id) const{

	lint num_of_ghosts = 0;

	for (int j = 0; j < _pDom->nZones; j++) {

		t_GhostLayer* glayer = _pGLayers[getPlainInd(zone_id, j)];
		num_of_ghosts += glayer->size();

	}

	return num_of_ghosts;

};

// calculate Offset in cell numbering for zone with id=zone_id_my
// so Cells[Offset] is the first ghost node 
// provided by zone with id zone_id_dnr
// NB: zone cell order: 
// [real cells][ghost from zone 0][ghost from zone 1]...[ghosts from zone N-1]
lint t_GhostMngBase::calcIndOffset(int zone_id_my, int zone_id_dnr) const{

	const t_Zone& Zne = _pDom->Zones[zone_id_my];

	// starting point of all ghosts
	lint offset = Zne.getnCellsReal();

	for (int j = 0; j < zone_id_dnr; j++) {

		const t_GhostLayer& glayer = *_pGLayers[j];
		offset += glayer.size();
	}

	return offset;

};

void t_GhostMngBase::getDonor(int iZone, int iCell, int& a_iZoneDnr, int& a_iCellDnr) const{

	// check that iCell is ghost
	const t_Zone& zne = _pDom->Zones[iZone];

	if (zne.isRealCell(iCell))
		hsLogError("t_GhostMngBase::getDonor: trying to find donor for real cell: iZone=%d, iCell=%d", 
			iZone, iCell);

	const t_GhostLayer* pLayer;

	for (int iZoneDnr = 0; iZoneDnr < _pDom->nZones; iZoneDnr++) {

		pLayer = _pGLayers[getPlainInd(iZone, iZoneDnr)];

		if (pLayer->size() > 0) {

			int idStart = calcIndOffset(iZone, iZoneDnr);
			int idEnd = idStart + pLayer->size() - 1;

			if (idStart <= iCell && iCell <= idEnd) {
				int shift = iCell - idStart;
				// donor 2 ghost data
				const t_Cell2GhostData& d2g = pLayer->data[shift];
				if (d2g.id_my != iCell) {
					hsLogError("t_GhostMngBase::getDonor: found entry, but cell 2 ghost data is wrong...");
					hsLogError("... tried to find donor for: iZone=%d, iCell=%d", iZone, iCell);
					hsLogError("... corresponding data in ghost layer: iCell=%d", d2g.id_my);
				}
				a_iZoneDnr = iZoneDnr;
				a_iCellDnr = d2g.id_dnr;
			}
		}
	}


};

void t_GhostMngBase::exchangeGeomData() {

	for (int i = 0; i < _pDom->nZones; i++) {

		for (int j = 0; j < _pDom->nZones; j++) {

			t_Zone& ZoneMy = _pDom->Zones[i];
			t_Zone& ZoneDnr = _pDom->Zones[j];

			t_GhostLayer* glayer = _pGLayers[getPlainInd(i, j)];

			for (int k = 0; k < glayer->size(); k++) {

				lint offset = calcIndOffset(i, j);
				ZoneMy.getCell(offset + k).Center = ZoneDnr.getCell(glayer->data[k].id_dnr).Center;

			}

		}
	}

}