#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"
#include "common_procs.h"

#include "bc_data.h"

#include "settings.h"

//
// Forward declarations
//
static bool parseConnectivity(t_CGNSContext& ctx);  // 1-to-1 connectivity

static bool check_BCs(t_CGNSContext& ctx);
static bool parseBCs(t_CGNSContext& ctx);     // boundary conditions
static bool parseVCs(t_CGNSContext& ctx);     // volume conditions (frozen zones)

void loadCells(t_CGNSContext& ctx);
static bool loadGridCoords(t_CGNSContext& ctx);

t_CellKind getElementKind(CG_ElementType_t cg_type) {

	t_CellKind cell_kind;

	if (cg_type == CG_TETRA_4) {
		cell_kind = t_CellKind::Tetra;
	}
	if (cg_type == CG_HEXA_8) {
		cell_kind = t_CellKind::Brick;
	}
	return cell_kind;

};

int read_cgns_mesh()
{

	cgsize_t isize[3];

	cgsize_t irmin, irmax, istart, iend;

	const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

	char zonename[CG_MAX_NAME_LENGTH], sectionname[CG_MAX_NAME_LENGTH];

	CGNS_ENUMT(ElementType_t) itype;

	int res = CG_OK;

	char szName[CG_MAX_NAME_LENGTH];  // names in CGNS file

	t_CGNSContext ctx;

	const char* gridFN = g_genOpts.strGridFN.c_str();

	if (cg_open(gridFN, CG_MODE_READ, &ctx.iFile) != CG_OK)
	{
		hsLogMessage("Can't open grid file '%s' for reading (%s)",
			gridFN, cg_get_error());
	}

	ctx.iBase = 1; // assume only one base

	//
	// Space dimensions
	//
	int dimCell = 0, dimPhys = 0;
	cg_base_read(ctx.iFile, ctx.iBase, szName, &dimCell, &dimPhys);
	if (dimCell != G_Domain.nDim)
	{
		hsLogMessage("The grid is not for %dD problems", G_Domain.nDim);
	}

	//
	// Number of zones
	//
	cg_nzones(ctx.iFile, ctx.iBase, &G_Domain.nZones);

	// assignZonesToProcs() will be here when multiblock is up
	// now we need only G_Domain.map_iZne2cgID

	G_Domain.map_iZne2cgID = new int[G_Domain.nZones];

	//
	// Default layout: one-to-one mapping of zones indices to CGNS zone IDs
	//
	for (int b = 0; b < G_Domain.nZones; ++b)
		G_Domain.map_iZne2cgID[b] = b + 1;

	// Allocate memory for whole computational domain
	G_Domain.Zones = new t_Zone[G_Domain.nZones];

	for (int i = 0; i < G_Domain.nZones; i++) G_Domain.Zones[i].setId(i);

	// Temporary zones data used on CGNS file parsing
	ctx.cgZones = new t_CGNSZone[G_Domain.nZones];


	for (int iZone = 0; iZone < G_Domain.nZones; iZone++) {

		const int& cgZneID = G_Domain.map_iZne2cgID[iZone];
		t_Zone& Zne = G_Domain.Zones[iZone];
		t_CGNSZone& cgZne = ctx.cgZones[iZone];

		hsLogMessage("Reading zone#%d", iZone);

		CG_ZoneType_t type;  cg_zone_type(ctx.iFile, ctx.iBase, cgZneID, &type);
		bool isOk = (type == CG_Unstructured) ? 1 : 0;
		if (!isOk) hsLogMessage("Only unstructured grids are supported");

		cg_zone_read(ctx.iFile, ctx.iBase, cgZneID, zonename, isize);

		// CGNS documentation: midlevel/structural.html#zone
		// isize = {NVertex, NCell3D, NBoundVertex}
		Zne.initialize(isize[0], isize[1]);

		const cgsize_t& nVerts = Zne.getnVerts();
		const cgsize_t& nCells = Zne.getnCells();

		hsLogMessage("Number of Verts:%d", nVerts);

		hsLogMessage("Number of Cells:%d", nCells);

		//hsLogMessage("Number of BC Verts:%d", isize[2]);

		ctx.map_ZneName2Idx[Zne.getName()] = iZone;
		
		// reading sections

		int nsections, index_sect, nbndry, iparent_flag;
		cgsize_t iparentdata;

		cg_nsections(ctx.iFile, ctx.iBase, cgZneID, &nsections);

		hsLogMessage("Number of sections:%d", nsections);

		for (index_sect = 1; index_sect <= nsections; index_sect++)
		{
			cgsize_t n_elems, n_verts_in_elem = 0;

			cg_section_read(ctx.iFile, ctx.iBase, cgZneID, index_sect, sectionname,
				&itype, &istart, &iend, &nbndry, &iparent_flag);
			printf("\nReading section info...\n");
			printf("   section name=%s\n", sectionname);
			printf("   section type=%s\n", ElementTypeName[itype]);
			printf("   istart,iend=%i, %i\n", (int)istart, (int)iend);

			// TODO: universal way to detect sections containing cells
			// for now detect by type of elements
			//if (itype != CG_HEXA_8 && itype != CG_TETRA_4) continue;

			if (itype == CG_HEXA_8) n_verts_in_elem = 8;
			if (itype == CG_TETRA_4) n_verts_in_elem = 4;

			n_elems = iend - istart + 1;

			// bc
			if (itype == CG_QUAD_4) n_verts_in_elem = 4;
			if (itype == CG_TRI_3) n_verts_in_elem = 3;

			// reading face patches
			if (itype == CG_QUAD_4 || itype == CG_TRI_3) {
				t_CGSection* pPatch = new t_CGSection(sectionname);

				t_FaceBCKind kind;
				if (G_BCList.getBCKindBySectName(sectionname, kind) == true)
					cgZne.addSection(pPatch, t_CGSectionKind::BC);
				else
					cgZne.addSection(pPatch, t_CGSectionKind::Abutted);

				pPatch->alloc(n_elems, n_verts_in_elem);

				hsLogMessage("   reading face patch : section %d, Type: %s\n",
					index_sect, ElementTypeName[itype]);

				cg_elements_read(ctx.iFile, ctx.iBase, cgZneID, index_sect, pPatch->get_buf_data(), \
					& iparentdata);

				pPatch->itype = itype;

				// debug output of section
				for (int i = 0; i < n_elems; i++) {
					for (int j = 0; j < n_verts_in_elem; j++)
						std::cout << pPatch->get_buf().get_val(i, j) << ";";
					std::cout << std::endl;
				}


				continue;
			}
			// reading elements
			if (itype == CG_HEXA_8 || itype == CG_TETRA_4) {

				t_CGSection* pNewCellSet = new t_CGSection();

				cgZne.addSection(pNewCellSet, t_CGSectionKind::Cell);

				pNewCellSet->alloc(n_elems, n_verts_in_elem);

				hsLogMessage("   reading elements : section %d, Type: %s\n",
					index_sect, ElementTypeName[itype]);

				cg_elements_read(ctx.iFile, ctx.iBase, cgZneID, index_sect, pNewCellSet->get_buf_data(), \
					& iparentdata);

				pNewCellSet->itype = itype;

				// debug output of section
				for (int i = 0; i < n_elems; i++) {
					for (int j = 0; j < n_verts_in_elem; j++)
						std::cout << pNewCellSet->get_buf().get_val(i, j) << ";";
					std::cout << std::endl;
				}
				continue;
			}

			hsLogMessage("Error: read_cgns_mesh(): unsupported section type");

		}

	}

	loadCells(ctx);

	// Read grid coordinates
	if (!loadGridCoords(ctx))
		return false;

	// Connectivity info
	if (!parseConnectivity(ctx))
		return false;

	G_Domain.makeVertexConnectivity();

	G_Domain.makeCellConnectivity();

	G_Domain.makeFaces();

	if (G_Domain.checkNormalOrientations()) 
		hsLogMessage("check Face Normal Orientations : Ok");
	else
		hsLogMessage("Error:checkNormalOrientations failed!");

	G_Domain.calcUnitOstrogradResid();

	// Volume conditions info (frozen zones)
	parseVCs(ctx);

	// checks with some cgns bc-specific funcs
	if (!check_BCs(ctx))
		return false;
	// update mesh with bc sets
	if (!parseBCs(ctx))
		return false;

	return 0;
}

void loadCells(t_CGNSContext& ctx) {

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		t_Zone& Zne = G_Domain.Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];
		
		cgsize_t NCellsCG = cgZne.countCells();
		if (NCellsCG != Zne.getnCells())
			hsLogError("loadCells: size mismatch of cells: %ld in CGNS Zone, %ld in Zone", 
			NCellsCG, Zne.getnCells());

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


static bool parse_1to1_connectivity_patch(const t_CGNSContext& ctx,
	const int iZne, const int cgPatchID, t_CGSection& cgFacePatch) {
	return true;
};

static bool parseConnectivity(t_CGNSContext& ctx) {

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		t_Zone& zne = G_Domain.Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		int n1to1 = 0;

		if (G_State.mpiRank == 0)
			cg_n1to1(ctx.iFile, ctx.iBase, cgZneID, &n1to1);

		int NPatchesAbut = 0;

		if (G_State.mpiRank == 0)
			cg_nconns(ctx.iFile, ctx.iBase, cgZneID, &NPatchesAbut);

		//ier = cg_nconns(int fn, int B, int Z, int* nconns)
		//cg_conn_info();
		//cg_conn_read();
		
		//MPI_Bcast(&cgZne.nPatches, 1, MPI_INT, 0, PETSC_COMM_WORLD);

		if (NPatchesAbut + n1to1 < 1 && G_Domain.nZones > 1) {
			hsLogError("Missing inter-zone 1-to-1 connectivity in '%s'#%d", zne.getName(), cgZneID);
			//return false;
		}

		if (n1to1 > 0 )
			hsLogError(
				"parseConnectivity: 1-to-1 connectivity detected, remake grid with generalized connectivity");

		if (NPatchesAbut != cgZne.getSectsAbut().size()) {
			hsLogError("parseConnectivity:Wrong number of abutted patches");
			//return false;
		}

		// Loop through connectivity patches
		for (int cgPatchId = 1; cgPatchId <= NPatchesAbut; ++cgPatchId)
		{
			t_CGSection& cgFacePatch = cgZne.getSectionAbut(cgPatchId - 1);

			char connectname[128];
			CG_GridLocation_t location;
			CG_GridConnectivityType_t connect_type;
			CG_PointSetType_t ptset_type;
			cgsize_t npnts;

			char donorname[128];
			CG_ZoneType_t donor_zonetype;
			CG_PointSetType_t donor_ptset_type;
			CG_DataType_t donor_datatype;
			cgsize_t ndata_donor[128];
			int ier = cg_conn_info(ctx.iFile, ctx.iBase, cgZneID, cgPatchId, connectname,
				&location, &connect_type,
				&ptset_type, &npnts, donorname,
				&donor_zonetype, &donor_ptset_type,
				&donor_datatype, ndata_donor);

			// ndata*3 for tris ndata*4 for quads
			cgsize_t pnts[64];
			cgsize_t donor_data[64];
			ier = cg_conn_read(ctx.iFile, ctx.iBase, cgZneID, cgPatchId,
				pnts, donor_datatype, donor_data);

			int bla = 1;

			//if (!parse_1to1_connectivity_patch(ctx, iZne, cgPatchId, cgFacePatch))
			//	return false;
		}

	}
	return true;
}

static bool check_BCs(t_CGNSContext& ctx) { 

	// bcs
	int nBCs;

	const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		t_Zone& Zne = G_Domain.Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		cg_nbocos(ctx.iFile, ctx.iBase, cgZneID, &nBCs);

		hsLogMessage("BC check .................Number of BC sets:%d", nBCs);

		for (int iBC = 1; iBC <= nBCs; ++iBC) {

			char szPatchName[CG_MAX_NAME_LENGTH];

			CG_BCType_t iBCtype;

			CG_PointSetType_t pntSetType;
			cgsize_t nPnts = 0; // number of points defining the BC region

								// Normals to the BC patch - UNUSED
			int iNorm[3]; // normal as index vector (computational coords)
			cgsize_t normListSize = 0;  CG_DataType_t normDataType; // normals in phys coords

			int nDatasets = 0; // number of datasets with additional info for the BC

			cg_boco_info(ctx.iFile, ctx.iBase, cgZneID, iBC,
				szPatchName, &iBCtype,
				&pntSetType, &nPnts,
				iNorm, &normListSize, &normDataType,
				&nDatasets
			);
			if (pntSetType != CG_PointRange && nPnts != 2)
			{
				hsLogMessage(
					"Boundary condition patch '%s'(#%d) of zone '%s'(#%d) isn't defined as point range",
					szPatchName, iBC, Zne.getName(), cgZneID);
				szPatchName[0] = 0x3;  // 'end of text' code -> error indicator
				return false;
			}
			cgsize_t idxRng[2];
			cg_boco_read(ctx.iFile, ctx.iBase, cgZneID, iBC, idxRng, nullptr);

			hsLogMessage("BC with name %s has idxrng_start=%d and idxrng_end=%d", szPatchName, idxRng[0], idxRng[1]);

			char szBC[CG_MAX_NAME_LENGTH] = "";

			cg_goto(ctx.iFile, ctx.iBase, "Zone_t", cgZneID, "ZoneBC", 0, "BC_t", iBC, NULL);
			if (cg_famname_read(szBC) == CG_OK)
			{
				// Read "Fam_Descr_Name" generated by Pointwise 16.03
				if (cg_goto(ctx.iFile, ctx.iBase, szBC, 0, NULL) == CG_OK)
				{
					char szDescrName[CG_MAX_NAME_LENGTH] = "";
					char* szDescrText = nullptr;
					if (cg_descriptor_read(1, szDescrName, &szDescrText) == CG_OK)
					{
						if (strcmp(szDescrName, "Fam_Descr_Name") == 0)
						{
							strcpy_s(szBC, szDescrText);
							cg_free(szDescrText);
						}
					}
				}
			}
			else
			{
				strcpy_s(szBC, szPatchName);
			}
		}
		hsLogMessage("BC check .................Ok");
		// ~bcs

	}

	return true; 
}
// parse BCs after they have already been read into ctx
// face lists in zones must be initialized too
static bool parseBCs(t_CGNSContext& ctx) {

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		t_Zone& Zne = G_Domain.Zones[iZne];
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

	return true;
};     // boundary conditions

static bool parseVCs(t_CGNSContext& ctx) { return true; }

static bool loadGridCoords(t_CGNSContext& ctx) {

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne) {

		double *x, *y, *z;

		t_Zone& Zne = G_Domain.Zones[iZne];

		const cgsize_t& nVerts = Zne.getnVerts();

		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		const char* coordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };

		x = new double[nVerts];
		y = new double[nVerts];
		z = new double[nVerts];

		cgsize_t irmin = 1; cgsize_t irmax = nVerts;

		cg_coord_read(ctx.iFile, ctx.iBase, cgZneID, "CoordinateX",
			CG_RealDouble, &irmin, &irmax, x);
		cg_coord_read(ctx.iFile, ctx.iBase, cgZneID, "CoordinateY",
			CG_RealDouble, &irmin, &irmax, y);
		cg_coord_read(ctx.iFile, ctx.iBase, cgZneID, "CoordinateZ",
			CG_RealDouble, &irmin, &irmax, z);

		for (int i = 0; i < Zne.getnVerts(); i++) {
			t_Vert& vert = Zne.getVert(i);
			vert.xyz[0] = x[i];
			vert.xyz[1] = y[i];
			vert.xyz[2] = z[i];
		}
	
		delete[] x, y, z;

	}

	return true;
}




