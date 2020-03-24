/*   Program read_grid_unst   */
/*
Reads simple 3-D unstructured grid from a CGNS file
(created using write_grid_unst.c).

The CGNS grid file 'grid_c.cgns' must already exist.

Example compilation for this program is (change paths if needed!):

cc -I ../.. -c read_grid_unst.c
cc -o read_grid_unst_c read_grid_unst.o -L ../../lib -lcgns

(../../lib is the location where the compiled
library libcgns.a is located)
*/

#include <stdio.h>
#include <iostream>

#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"
#include "common_procs.h"

//
// Forward declarations
//
static bool parseConnectivity(t_CGNSContext& ctx);  // 1-to-1 connectivity
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

	// READ X, Y, Z GRID POINTS FROM CGNS FILE
	// open CGNS file for read-only
	int res = CG_OK;

	char szName[CG_MAX_NAME_LENGTH];  // names in CGNS file

	t_CGNSContext ctx;

	char gridFN[] = "test_case/box-hexa-simple-2blk.cgns";
	//char gridFN[] = "test_case/box-tetra-simple.cgns";
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

		// isize[0] is nVerts, isize[1] is nCells
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
			if (itype != CG_HEXA_8 && itype != CG_TETRA_4) continue;

			if (itype == CG_HEXA_8) n_verts_in_elem = 8;
			if (itype == CG_TETRA_4) n_verts_in_elem = 4;

			n_elems = iend - istart + 1;

			// bc
			//if (itype == CG_QUAD_4) n_verts_in_elem = 4;
			//if (itype == CG_TRI_3) n_verts_in_elem = 3;

			cgZne.cells.allocate(n_elems, n_verts_in_elem);

			hsLogMessage("   reading element data for %s\n", ElementTypeName[itype]);
			cg_elements_read(ctx.iFile, ctx.iBase, cgZneID, index_sect, cgZne.cells.data(), \
				&iparentdata);

			cgZne.itype = itype;

			// debug output of sections
			for (int i = 0; i < n_elems; i++) {
				for (int j = 0; j < n_verts_in_elem; j++)
					std::cout << cgZne.cells.get_val(i, j) << ";";
				std::cout << std::endl;
			}

		}

	}

	loadCells(ctx);

	// Read grid coordinates
	if (!loadGridCoords(ctx))
		return false;

	G_Domain.makeVertexConnectivity();

	G_Domain.makeCellConnectivity();

	G_Domain.makeFaces();

	// Volume conditions info (frozen zones)
	parseVCs(ctx);

	// Connectivity info
	// Updates {zne,cgZne}.{is,ie,js,je,ks,ke}, zne.{nx,ny,nz}
	if (!parseConnectivity(ctx))
		return false;

	// Boundary conditions info
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

		if (cgZne.cells.nRows != Zne.getnCells()) hsLogMessage("loadCells: size mismatch of CGNS and working Zones");

		for (int i = 0; i < Zne.getnCells(); i++) {

			t_Cell& cell = Zne.getCell(i);

			cell.setKind(getElementKind(cgZne.itype));

			// additional check
			if (cell.NVerts != cgZne.cells.nCols) hsLogMessage("loadCells: wrong number of vertices in cell");
			// cgns IDs are 1-based, we use zero-based ids
			for (int j = 0; j < cell.NVerts; j++) {
				lint Vert_ID = cgZne.cells.get_val(i, j) - 1;
				cell.pVerts[j] = &(Zne.getVert(Vert_ID));
			}
		}
		std::cout << "______________________Debug, Zone Verts:\n";
		// debug output of sections
		for (int i = 0; i < Zne.getnCells(); i++) {
			const t_Cell& cell = Zne.getCell(i);
			for (int j = 0; j < cell.NVerts; j++)
				std::cout << Zne.getCell(i).getVert(j).Id<< ";";
			std::cout << std::endl;
		}

	}

};



static bool parseConnectivity(t_CGNSContext& ctx) { return true; }

static bool parseBCs(t_CGNSContext& ctx) { 

	// bcs
	int nBCs;

	const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

	for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
		t_Zone& Zne = G_Domain.Zones[iZne];
		t_CGNSZone& cgZne = ctx.cgZones[iZne];

		cg_nbocos(ctx.iFile, ctx.iBase, cgZneID, &nBCs);

		hsLogMessage("Number of BC sets:%d", nBCs);

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

		// ~bcs

	}

	
	
	return true; }

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
			vert.xyz.x = x[i];
			vert.xyz.y = y[i];
			vert.xyz.z = z[i];

		}
	
		delete[] x, y, z;

	}

	return true;
}




