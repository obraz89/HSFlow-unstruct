#include "common_data.h"

#include "flow_common.h"


// TODO: remove this dependency, how to store flux for different models?
#include "flow_model.h"

#include "ghost_manager.h"

#include <fstream>

#include "CGNS-ctx.h"

#include "bc_data.h"

t_CellKind getElementKind(CG_ElementType_t cg_type) {

	t_CellKind cell_kind = t_CellKind::None;

	if (cg_type == CG_TETRA_4) {
		cell_kind = t_CellKind::Tetra;
	}
	if (cg_type == CG_HEXA_8) {
		cell_kind = t_CellKind::Brick;
	}
	return cell_kind;

};

void t_Domain::initializeFromCtx() {

	nZones = G_CGNSCtx.nZones;

	// assignZonesToProcs() will be here when multiblock is up
// now we need only G_Domain.map_iZne2cgID

	map_iZne2cgID = new int[G_Domain.nZones];

//
// Default layout: one-to-one mapping of zones indices to CGNS zone IDs
//
	for (int b = 0; b < nZones; ++b)
		map_iZne2cgID[b] = b + 1;

	Zones = new t_Zone[nZones];

	for (int i = 0; i < nZones; i++) { 

		t_Zone& zne = Zones[i];
		int cg_id = map_iZne2cgID[i];
		const t_CGNSZone& cgZne = G_CGNSCtx.cgZones[i];
		
		zne.setIdGlob(i); 

		zne.setName(cgZne.getName());


	
	}

	// bunch of code from read_cgns_mesh()
	{


		// get sizes of zones and read real cells from cgns ctx
		loadCells();

		// set up connections from verts to real cells
		makeVertexConnectivity();

		G_GhostManager.setDom(*this);
		G_GhostManager.initialize(G_CGNSCtx);

		makeCellConnectivity();

		makeFaces();


		// update mesh with bc sets
		loadBCs();

		if (checkNormalOrientations())
			hsLogMessage("check Face Normal Orientations : Ok");
		else
			hsLogMessage("Error:checkNormalOrientations failed!");

		calcUnitOstrogradResid();

		// Volume conditions info (frozen zones)
		//parseVCs(ctx);

	}



}

void t_Domain::loadCells() {

	const t_CGNSContext& ctx = G_CGNSCtx;

	for (int iZne = 0; iZne < nZones; ++iZne)
	{
		const int& cgZneID = map_iZne2cgID[iZne];
		t_Zone& Zne = Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		cgsize_t NCellsCG = cgZne.countCells();

		// check that ncells is ok
		if (NCellsCG != cgZne.getNCells())
			hsLogError("loadCells: number of cells in cgZne is different from what was read from cgns file!");

		cgsize_t NGhosts = ctx.getNumOfGhostsForZone(cgZneID);
		// Zone stores real cells + some ghost cells

		cgsize_t NCellsTot = NCellsCG + NGhosts;

		Zne.initialize(cgZne.getNVerts(), NCellsCG, NCellsTot);

		if (NCellsCG > Zne.getnCellsTot())
			hsLogError("loadCells: failed to initialize Zne: %ld cells in CGNS Zone, %ld cells in Zone",
				NCellsCG, Zne.getnCellsTot());

		// filling in real cells
		int iCell = 0;

		for (int i = 0; i < cgZne.getSectsCell().size(); i++) {

			const t_CGSection& cg_cells = cgZne.getSectionCell(i);

			for (int j = 0; j < cg_cells.get_buf().nRows; j++) {

				t_Cell& cell = Zne.getCell(iCell);

				cell.setKind(getElementKind(cg_cells.itype));

				for (int k = 0; k < cg_cells.get_buf().nCols; k++) {
					// cgns IDs are 1-based, we use zero-based ids
					lint Vert_ID = cg_cells.get_buf().get_val(j, k) - 1;
					cell.pVerts[k] = &(Zne.getVert(Vert_ID));
				}

				iCell++;

			}
		}
		// load grid coords
		for (int i = 0; i < cgZne.getNVerts(); i++) {
			t_Vert& vert = Zne.getVert(i);
			vert.xyz[0] = cgZne.getXCoords()[i];
			vert.xyz[1] = cgZne.getYCoords()[i];
			vert.xyz[2] = cgZne.getZCoords()[i];
		}

		//std::cout << "______________________Debug, Zone Verts:\n";
		// debug output of sections
		//for (int i = 0; i < Zne.getnCells(); i++) {
		//	const t_Cell& cell = Zne.getCell(i);
		//	for (int j = 0; j < cell.NVerts; j++)
		//		std::cout << cell.getVert(j).Id<< ";";
		//	std::cout << std::endl;
		//}

	}

};

// parse BCs after they have already been read into ctx
// face lists in zones must be initialized too
void t_Domain::loadBCs() {

	const t_CGNSContext& ctx = G_CGNSCtx;

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int cgZneID = iZne + 1;
		t_Zone& Zne = Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		for (int ipatch = 0; ipatch < cgZne.getSectsBC().size(); ipatch++) {

			const t_CGSection& fpatch_cg = cgZne.getSectionBC(ipatch);

			t_FaceBCKind bc_kind;
			bool ok = G_BCList.getBCKindBySectName(fpatch_cg.name, bc_kind);
			//ok if the patch is a bc patch (not a zone-2-zone patch)
			if (ok) {
				int nF = fpatch_cg.get_buf().nRows;
				int NVertsInFace = fpatch_cg.get_buf().nCols;
				t_Face* flist = new t_Face[nF];

				for (int iF = 0; iF < nF; iF++) {

					t_Face& face = flist[iF];

					face.BCKind = bc_kind;

					face.NVerts = NVertsInFace;

					for (int j = 0; j < NVertsInFace; j++) {

						cgsize_t iVert = fpatch_cg.get_buf().get_val(iF, j);
						// cg Id is 1-based, we make our 0-based
						face.pVerts[j] = Zne.getpVert(iVert - 1);
					}

				}

				Zne.updateFacesWithBCPatch(flist, nF);

				delete[] flist;
			}
		}
	}
};     // boundary conditions

void t_Domain::initializeFlow() {

	t_ConsVars cvs = calcConsVarsInf();

	// set real cell values
	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		for (int i = 0; i < zne.getnCellsReal(); i++) {

			pcell = zne.getpCell(i);
			pcell->ConsVars = cvs;
		}

	}
	// set ghost values
	G_GhostManager.exchangeCSV();

}

void t_Domain::dump_flow() {

	std::string fn("dump_flow.txt");

	std::ofstream ofstr(fn);

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		ofstr << "=========Zone #" << iZone << "===========\n";

		for (int i = 0; i < zne.getnCellsTot(); i++) {

			pcell = zne.getpCell(i);
			ofstr <<"cell #"<<i << pcell->ConsVars.to_str();
		}

	}

	ofstr.flush();

}

void t_Domain::dump_geom() {



}