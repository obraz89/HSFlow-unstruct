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
/* cgnslib.h file must be located in directory specified by -I during compile: */
#include "cgnslib.h"

#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"

int read_cgns_mesh()
{
	double *x, *y, *z;
	//float x[21 * 17 * 9], y[21 * 17 * 9], z[21 * 17 * 9];

	cgsize_t *isize, *ielem;
	//cgsize_t isize[3][1], ielem[20 * 16 * 8][8];
	int index_file, index_base, index_zone;
	cgsize_t irmin, irmax, istart, iend;
	int nsections, index_sect, nbndry, iparent_flag;
	cgsize_t iparentdata;
	char zonename[33], sectionname[33];
	CGNS_ENUMT(ElementType_t) itype;

	// READ X, Y, Z GRID POINTS FROM CGNS FILE
	// open CGNS file for read-only
	int res = CG_OK;

	const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

	char szName[CG_MAX_NAME_LENGTH];  // names in CGNS file

	t_CGNSContext ctx;

	char gridFN[] = "test_case/box_10x10x10.cgns";
	if (cg_open(gridFN, CG_MODE_READ, &ctx.fileID) != CG_OK)
	{
		hsLogMessage("Can't open grid file '%s' for reading (%s)",
			gridFN, cg_get_error());
	}

	ctx.iBase = 1; // assume only one base

	//
	// Space dimensions
	//
	int dimCell = 0, dimPhys = 0;
	cg_base_read(ctx.fileID, ctx.iBase, szName, &dimCell, &dimPhys);
	if (dimCell != G_Domain.nDim)
	{
		hsLogMessage("The grid is not for %dD problems", G_Domain.nDim);
	}

	//
	// Number of zones
	//
	cg_nzones(ctx.fileID, ctx.iBase, &G_Domain.nZones);
	if (G_Domain.nZones != 1)
	{
		hsLogMessage("The domain contains not a single zone (%d)", G_Domain.nZones);
	}
	// Allocate memory for whole computational domain
	G_Domain.Zones = new t_Zone[G_Domain.nZones];

	// Temporary zones data used on CGNS file parsing
	ctx.cgZones = new t_CGNSZone[G_Domain.nZones];

	const int cgZneID = 1;

	CG_ZoneType_t type;  cg_zone_type(ctx.fileID, ctx.iBase, cgZneID, &type);
	bool isOk = (type == CG_Unstructured) ? 1 : 0;
	if (!isOk) hsLogMessage("Only unstructured grids are supported");


	/*
	// get zone size (and name - although not needed here)
	cg_zone_read(index_file, index_base, index_zone, zonename, isize[0]);
	// lower range index
	irmin = 1;
	// upper range index of vertices
	irmax = isize[0][0];
	// read grid coordinates
	cg_coord_read(index_file, index_base, index_zone, "CoordinateX",
		CGNS_ENUMV(RealSingle), &irmin, &irmax, x);
	cg_coord_read(index_file, index_base, index_zone, "CoordinateY",
		CGNS_ENUMV(RealSingle), &irmin, &irmax, y);
	cg_coord_read(index_file, index_base, index_zone, "CoordinateZ",
		CGNS_ENUMV(RealSingle), &irmin, &irmax, z);
	// find out how many sections
	cg_nsections(index_file, index_base, index_zone, &nsections);
	printf("\nnumber of sections=%i\n", nsections);
	// read element connectivity
	for (index_sect = 1; index_sect <= nsections; index_sect++)
	{
		cg_section_read(index_file, index_base, index_zone, index_sect, sectionname,
			&itype, &istart, &iend, &nbndry, &iparent_flag);
		printf("\nReading section data...\n");
		printf("   section name=%s\n", sectionname);
		printf("   section type=%s\n", ElementTypeName[itype]);
		printf("   istart,iend=%i, %i\n", (int)istart, (int)iend);
		if (itype == CGNS_ENUMV(HEXA_8))
		{
			printf("   reading element data for this element\n");
			cg_elements_read(index_file, index_base, index_zone, index_sect, ielem[0], \
				&iparentdata);
		}

		if (itype == CG_TETRA_4) {
			printf("   reading element data for tetra\n");
			cg_elements_read(index_file, index_base, index_zone, index_sect, ielem[0], \
				&iparentdata);
		}
		else
		{
			printf("   not reading element data for this element\n");
		}
	}
	// close CGNS file
	cg_close(index_file);
	printf("\nSuccessfully read unstructured grid from file grid_c.cgns\n");
	printf("   for example, element 1 is made up of nodes: %i, %i, %i, %i, %i, %i, %i, %i\n",
		(int)ielem[0][0], (int)ielem[0][1], (int)ielem[0][2], (int)ielem[0][3],
		(int)ielem[0][4], (int)ielem[0][5], (int)ielem[0][6], (int)ielem[0][7]);
	printf("   x,y,z of node 357 are: %f, %f, %f\n", x[357], y[357], z[357]);
	printf("   x,y,z of node 1357 are: %f, %f, %f\n", x[1357], y[1357], z[1357]);
	return 0;
	*/
	return 0;
}
