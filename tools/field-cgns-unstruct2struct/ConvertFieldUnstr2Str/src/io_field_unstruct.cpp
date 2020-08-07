/**********************************
// Name:        io-field.cpp
// Purpose:     Field input-output routines, cgns format
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
**********************************/

#include "mpi.h"

#include "cgnslib.h"
// TODO: make cgns wrapper as in HSFlow
#if CG_BUILD_SCOPE == 0
#undef  CG_BUILD_SCOPE
#define CG_BUILD_SCOPE 1
#endif

#include "common_data_unstruct.h"

#include "logging.h"

#include "settings.h"

#include "io_field_unstruct.h"

#include "dom_base_unstruct.h"

#include "CGNS-ctx-unst.h"

#include <cassert>
//-----------------------------------------------------------------------------

const char g_szCGBase[] = "HSFlow-unstruct";

//-----------------------------------------------------------------------------

//
// Forward declarations
//
//-----------------------------------------------------------------------------

//
// Helper classes
//



/**
 *  Read zone data from the opened CGNS file
 *
 *  @param[in] fileID - ID of the opened CGNS file
 *  @param[in] iBase - base number in the file
 *  @param[in] idxZne - 0-based internal zone index
 *  @param[in] nxyz - dimensions of the zone, without ghosts!!!
 *
 *  @param[out] grid - node coordinates of the loaded block's grid
 *  @param[out] field - loaded field data
 *
 *  @return `true` if succeeded and `false` otherwise
**/
bool read_zone_cgns(const int fileID, const int iBase, const int idxZne,
					TpakArraysDyn<double>& field)
{
	const int& cgZneID = g_pDomUnst->map_iZne2cgID[idxZne];
	const t_Zone& Zne = g_pDomUnst->Zones[idxZne];

	// Get zone size and name
	char szZone[33];
	cgsize_t isize[3];	// nVerts, nCells, nBoundVerts

	if (cg_zone_read(fileID, iBase, cgZneID, szZone, isize) != CG_OK) {
		hsLogError("Can't read CGNS zone #%d ( %s )", cgZneID, cg_get_error());
		return false;
	}

	if (strcmp(szZone, Zne.getName()) != 0) {
		hsLogWarning("Inconsistent zone #%d names: '%s' -> '%s'",
			cgZneID, szZone, Zne.getName());
	}

	// Block size without ghosts
	if (isize[1] != Zne.getnCellsReal())
	{
		hsLogError("Inconsistent zone '%s'#%d dimensions (ok if this is Fluent field): %d <-> %d",
			szZone, cgZneID,
			isize[1], Zne.getnCellsReal());
		//return false;
	}
	// Indexes faces
	cgsize_t irmin = 1;
	cgsize_t irmax = Zne.getnVerts();

	//
	// Get solution info
	// FIXME: flow assumed existing
	//
	int iSol = 1;
	{
		CG_GridLocation_t loc;   char cgName[33];
		if (cg_sol_info(fileID, iBase, cgZneID, iSol, cgName, &loc) != CG_OK)
		{
			hsLogError("Can't read flow info from zone '%s'#%d ( %s )",
				szZone, cgZneID, cg_get_error());
			return false;
		}

		if (loc != CG_CellCenter)
		{
			hsLogError("CGNS: GridLocation must be CellCenter");
			return false;
		}
	}

	//
	// Read functions
	//
	assert(field.size() > 0);
	std::vector<std::string> funcNames = g_pDomUnst->getFuncNamesIO();

	irmin = 1;
	irmax = Zne.getnCellsReal();

	for (int idf = 0; idf < funcNames.size(); idf++) {

		double* __restrict U = field[idf];

		const char* fun_name = funcNames[idf].c_str();

		int r = cg_field_read(fileID, iBase, cgZneID, iSol, fun_name,
			CG_RealDouble, &irmin, &irmax, U);

		if (r != CG_OK && r != CG_NODE_NOT_FOUND) {
			hsLogError("Can't read '%s' from zone '%s'#%d ( %s )",
				fun_name, szZone, cgZneID, cg_get_error());
			return false;
		}

	}

	return true;
}
