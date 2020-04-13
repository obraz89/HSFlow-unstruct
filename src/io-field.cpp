/**********************************
// Name:        io-field.cpp
// Purpose:     Field input-output routines, cgns format
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
**********************************/



#include "cgnslib.h"
// TODO: make cgns wrapper as in HSFlow
#if CG_BUILD_SCOPE == 0
#undef  CG_BUILD_SCOPE
#define CG_BUILD_SCOPE 1
#endif

//#include "mpi.h"

#include "common_data.h"

#include "logging.h"
#include "settings.h"

#include "io-field.h"
#include "mesh.h"

#include "flow_params.h"

#include "CGNS-ctx.h"

#include <cassert>
//-----------------------------------------------------------------------------

const char g_szCGBase[] = "HSFlow-unstruct";

//-----------------------------------------------------------------------------

//
// Forward declarations
//
bool writeMetaInfoToCGNS(const int fileID, const int iBase, const short time_layer);
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
	const int& cgZneID = G_pMesh->map_iZne2cgID[idxZne];
	const t_Zone& Zne = G_pMesh->Zones[idxZne];

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
	if (isize[1]!=Zne.getnCellsReal())
	{
		hsLogError("Inconsistent zone '%s'#%d dimensions: %d <-> %d",
			szZone, cgZneID,
			isize[1], Zne.getnCellsReal());
		return false;
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
	std::vector<std::string> funcNames = G_pMesh->getFuncNamesIO();

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


/**
 *  Saves field data to file
 *
 * @param[in] fileName - file name to save field data
 * @param[in] gridFileName - file name for grid, if empty embed grid into field file
 * @param[in] time_layer - time layer to save field from
 *                         0, 1 or 2 -> current, previous, pre-previous
 * @param[in] isDouble - use double or single precision for field arrays
 *
 * @return `true` if succeeded and `false` otherwise
**/
bool saveField(const std::string& fileName, const std::string& gridFileName,
	const short time_layer)
{

	int f = -1, fGrid = -1;  // file descriptors for field & grid
	std::string fn_grid;  // grid file name relative to the field
	int iBase = -1, iBaseGrid = -1;  // CGNS base ids

	TLogSyncGuard logGuard;
	short ok = 1;  // error code for MPI broadcasting
	if (G_State.mpiRank == 0) do
	{
		ok = 0;

		const std::string path_field = g_CASE_RESULTS_DIR + fileName + ".cgns";

		if (cg_open(path_field.c_str(), CG_MODE_WRITE, &f) != CG_OK)
		{
			hsLogError("Can't open file '%s' for writing ( %s )",
				path_field.c_str(), cg_get_error());
			break;
		}
		hsLogMessage("\n* Saving field to file '%s'...", path_field.c_str());

		fGrid = f;
		std::string path_grid = path_field;

		if (!gridFileName.empty())
		{
			//
			// Create separate grid file
			//
			fn_grid = gridFileName + ".cgns";
			path_grid = g_CASE_RESULTS_DIR + fn_grid;

			static bool needToSave = true;
			if (needToSave)
			{
				if (cg_open(path_grid.c_str(), CG_MODE_WRITE, &fGrid) == CG_OK)
					needToSave = false;
				else {
					hsLogError("Can't open file '%s' for writing ( %s )",
						path_grid.c_str(), cg_get_error());
					// this often happens under Windows, when HDF-lib erroneously
					// doesn't release grid file after initial field loading

					// Don't stop & proceed to field saving,
					// the grid may be obtained offline
					fGrid = -1;
				}
			}
			else
				fGrid = -1;
		}

		// Create base
		if (cg_base_write(f, g_szCGBase, G_pMesh->nDim, G_pMesh->nDim, &iBase) != CG_OK)
		{
			hsLogError("Can't write base into '%s' ( %s )",
				path_field.c_str(), cg_get_error());
			break;
		}

		// Write Reference state, Equations info, time, etc
		if (!writeMetaInfoToCGNS(f, iBase, time_layer))
			break;

		iBaseGrid = iBase;
		if (fGrid != f && fGrid >= 0) // grid in separate file
		{
			if (cg_base_write(fGrid, g_szCGBase, G_pMesh->nDim, G_pMesh->nDim, &iBaseGrid) != CG_OK)
			{
				hsLogError("Can't write CGNS base node into '%s' ( %s )",
					path_grid.c_str(), cg_get_error());
				break;
			}
		}

		ok = 1;
	} while (false); //if( G_State.mpiRank == 0 )

	//MPI_Bcast(&ok, 1, MPI_SHORT, 0/*root*/, PETSC_COMM_WORLD);
	if (!ok)
		return false;

	// Inform all ranks if we need grid coordinates or not
	//MPI_Bcast(&fGrid, 1, MPI_INT, 0/*root*/, PETSC_COMM_WORLD);

	//
	// Sorted 1-based CGNS zone IDs and corresponding internal 0-based zone indices
	static int* map_cgID2iZne = nullptr;
	if (!map_cgID2iZne)
	{
		map_cgID2iZne = new int[G_pMesh->nZones];
		for (int b = 0; b < G_pMesh->nZones; ++b)
			map_cgID2iZne[G_pMesh->map_iZne2cgID[b] - 1] = b;
	}

	//
	// Write field & grid data
	//
	for (int cgZneID = 1; cgZneID <= G_pMesh->nZones; ++cgZneID)
	{
		const int& zi = map_cgID2iZne[cgZneID - 1];
		int iZone = -1;   int iZoneGrid = -1;

		t_Zone& zne = G_pMesh->Zones[zi];

		if (G_State.mpiRank == 0)
		{
			// Zone size packed in CGNS format
			cgsize_t isize[3];

			// assign isize here

			isize[0] = zne.getnVerts();
			isize[1] = zne.getnCellsReal();
			isize[2] = 0;

			// Zone for fields
			cg_zone_write(f, iBase, zne.getName(), isize, CG_Unstructured, &iZone);

			// Zone for grid coords
			iZoneGrid = iZone;
			if (fGrid != f && fGrid >= 0) // grid in separate file
			{
				cg_zone_write(fGrid, iBaseGrid, zne.getName(), isize, CG_Unstructured, &iZoneGrid);
			}

			//
			// Make link to the grid in separate file
			//
			if (fGrid != f)
			{
				// CGNS-node path inside the separate grid file
				const char* label = "GridCoordinates";
				const std::string& path = hs_string_format("/%s/%s/%s", g_szCGBase, zne.getName(), label);

				// Move CGNS file position to the current zone then write link
				cg_goto(f, iBase, "Zone_t", iZone, NULL);
				cg_link_write(label, fn_grid.c_str(), path.c_str());
			}
		}


		t_ArrDbl Vals;
		Vals.alloc(zne.getnVerts());

		//
		// Grid coords
		//
		if (fGrid >= 0)
		{
			for (int iCoord = 0; iCoord < 3; iCoord++) {
				const char* name = g_cgCoordNames[iCoord];
				int iCoordCG;
				G_pMesh->getDataAsArr(name, zi, Vals);
				if (cg_coord_write(fGrid, iBaseGrid, iZoneGrid, 
					CG_RealDouble, name, Vals.data(), &iCoordCG) != CG_OK) {
						hsLogError("Can't write grid coords in zone %s#%d ( %s )",
							zne.getName(), iZone, cg_get_error());
				};
			}
		} // if( fGrid >=0 )

		//
		// TODO: write bcs here for post-processing
		//

		// 2) write cells here
		const t_CGNSZone& cgZne = G_CGNSCtx.cgZones[zi];
		const std::vector<t_CGSection*> SectsCell = cgZne.getSectsCell();
		for (int iSect = 0; iSect < SectsCell.size(); iSect++) {

			const t_CGSection& Sect = *SectsCell[iSect];
			int sect_ind_cg_output = -1;
			const cgsize_t* data = Sect.get_buf_data();

			if (cg_section_write(f, iBase, iZone,
				Sect.name.c_str(), Sect.itype, Sect.id_start, Sect.id_end,
				0, data, &sect_ind_cg_output) != CG_OK) {
				hsLogError("Can't write cell section %s in zone %s#%d ( %s )",
					Sect.name.c_str(), zne.getName(), iZone, cg_get_error());
			};

		}

		// 3) write flow solution

		std::vector<std::string> flow_vars = G_pMesh->getFuncNamesIO();

		int iSol;

		cg_sol_write(f, iBase, iZone, "FlowSolution", CG_CellCenter, &iSol);


		t_ArrDbl flow_sol_data;
		flow_sol_data.alloc(zne.getnCellsReal());

		for (int iFlow = 0; iFlow < flow_vars.size(); iFlow++) {

			const char* flow_sol_name = flow_vars[iFlow].c_str();
			G_pMesh->getDataAsArr(flow_sol_name, zi, flow_sol_data);

			int iFlowCG;

			if (cg_field_write(f, iBase, iZone, iSol,
				CG_RealDouble, flow_sol_name,
				flow_sol_data.data(), &iFlowCG) != CG_OK) {

				hsLogError("Can't write %s in zone %s#%d ( %s )",
					flow_sol_name, zne.getName(), iZone, cg_get_error());
			};

		}


		//
		// Field data
		//
		// example code from structured=======================================================
		/*
		int iSol = -1;
		if (G_State.mpiRank == 0)
			cg_sol_write(f, iBase, iZone, "FlowSolution", CG_Vertex, &iSol);

		const double* srcU = (!zne.isFrozen) ? zne.UU[time_layer] : zne.U;
		for (int fun = 0; fun < G_Domain.nu; ++fun)
		{
			const int mpiTag = 's' + 'o' + 'l' + fun;

			if (G_Domain.bs <= zi && zi <= G_Domain.be)
			{
				//
				// Worker mpi-rank -> fill-in and send
				//
				size_t c = 0;
				for (int k = ks1; k <= ke1; ++k) {
					for (int j = js1; j <= je1; ++j) {
						for (int i = is1; i <= ie1; ++i)
						{
							int pos = G_Domain.nu * (zne.flatIdx(i, j, k) - 1) + fun;
							Vals.set(c++, srcU[pos]);
						}
					}
				}

				if (G_State.mpiRank != 0)
				{
					// NB: Don't use MPI_Send - it may flood the root MPI rank such that MPI_Recv fails
					MPI_Ssend(Vals.data(), nxyz0, Vals.type().mpi, 0, mpiTag, PETSC_COMM_WORLD);
				}
			}

			if (G_State.mpiRank == 0)
			{
				//
				// Root mpi-rank -> receive and write
				//
				const int& rankSrc = G_State.map_zone2rank[zi];
				if (rankSrc != 0) // don't receive from myself
				{
					// if writing was failed previously, then do a dummy recieve to unblock sender
					MPI_Recv(Vals.data(), (ok ? nxyz0 : 0), Vals.type().mpi,
						rankSrc, mpiTag, PETSC_COMM_WORLD,
						MPI_STATUS_IGNORE  // don't use NULL as status, MPI_Ssend may get stuck (i.e. in Intel MPI 5.0.1)
					);
				}

				int iField = -1;
				const char* name = G_Domain.phys->vecFuncNames[fun].c_str();
				if (ok)
					if (cg_field_write(f, iBase, iZone, iSol, Vals.type().cgns, name, Vals.data(), &iField) != CG_OK)
					{
						hsLogError("Can't write %s in zone %s#%d ( %s )",
							name, zne.szName, iZone, cg_get_error());
						ok = 0;  // don't return and continue recieving data from other ranks
					}
			}
		} // for( fun )
		*/
		// example code from structured=======================================================

	}  // Loop through zones

	if (G_State.mpiRank == 0)
	{
		cg_close(f);

		if (fGrid != f && fGrid >= 0)
			cg_close(fGrid);
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	//hsLogWTime();

	return ok;
}


/**
 *  Saves meta information about the field to the openned CGNS file
 *
 * @param[in] f     - ID of the opened CGNS file
 * @param[in] iBase - base in the CGNS file (mostly = 1)
 * @param[in] time_layer - time layer to save field from
 *                         0, 1 or 2 -> current, previous, pre-previous
 * @return `true` if succeded and `false` otherwise
**/
static bool writeMetaInfoToCGNS(const int f, const int iBase, const short time_layer)
{
	// Indicate class of the solution data
	cg_goto(f, iBase, NULL);
	cg_dataclass_write(CG_NormalizedByUnknownDimensional);

	//
	// Reference state
	//
	cg_goto(f, iBase, NULL);
	cg_state_write("ReferenceQuantities");
	{
		const cgsize_t one = 1;
		int count = 0;  // reference value number

		// Mach
		double Mach = G_FreeStreamParams.getMach();
		if (!isnan(Mach))
		{
			cg_goto(f, iBase, "ReferenceState_t", 1, NULL);
			cg_array_write("Mach", CG_RealDouble, 1, &one, &Mach);
			cg_goto(f, iBase, "ReferenceState_t", 1, "DataArray_t", ++count, NULL);
			cg_dataclass_write(CG_NondimensionalParameter);
		}

		// Re
		/*
		if (!isnan(refQ.Re))
		{
			cg_goto(f, iBase, "ReferenceState_t", 1, NULL);
			cg_array_write("Reynolds", CG_RealDouble, 1, &one, &refQ.Re);
			cg_goto(f, iBase, "ReferenceState_t", 1, "DataArray_t", ++count, NULL);
			cg_dataclass_write(CG_NondimensionalParameter);
		}
		*/

		// Temperature
		double TinfDim = G_FreeStreamParams.getTinfDim();
		if (!isnan(TinfDim))
		{
			cg_goto(f, iBase, "ReferenceState_t", 1, NULL);
			cg_array_write("Temperature", CG_RealDouble, 1, &one, &TinfDim);
			cg_goto(f, iBase, "ReferenceState_t", 1, "DataArray_t", ++count, NULL);
			cg_dataclass_write(CG_Dimensional);
		}
	}

	//
	// Flow equation set
	//
	cg_goto(f, iBase, NULL);
	/*
	if (cg_equationset_write(G_pMesh->nDim) == CG_OK)
	{
		const cgsize_t one = 1;
		std::map<std::string, double>::const_iterator iterPrm;

		// Gas model
		cg_goto(f, iBase, "FlowEquationSet_t", 1, NULL);
		cg_model_write("GasModel_t", CG_CaloricallyPerfect);

		// TODO: Chemical kinetics
		iterPrm = G_Domain.mapCasePrms_real.find("gamma");
		if (iterPrm != G_Domain.mapCasePrms_real.end())
		{
			// Cp/Cv
			const double gamma = iterPrm->second;
			cg_goto(f, iBase, "FlowEquationSet_t", 1, "GasModel_t", 1, NULL);
			cg_array_write("SpecificHeatRatio", CG_RealDouble, 1, &one, &gamma);
			cg_goto(f, iBase, "FlowEquationSet_t", 1, "GasModel_t", 1, "DataArray_t", 1, NULL);
			cg_dataclass_write(CG_NondimensionalParameter);

			// Ideal Gas Constant R = 1/(gamma*M^2) in nondimensional case
			{
				const double& M = G_Domain.phys->ref_state.M;
				const double R = 1. / (gamma * M * M);

				cg_goto(f, iBase, "FlowEquationSet_t", 1, "GasModel_t", 1, NULL);
				cg_array_write("IdealGasConstant", CG_RealDouble, 1, &one, &R);
				cg_goto(f, iBase, "FlowEquationSet_t", 1, "GasModel_t", 1, "DataArray_t", 2, NULL);
				cg_dataclass_write(CG_NondimensionalParameter);
			}
		}


		// Thermal Conductivity Model
		cg_goto(f, iBase, "FlowEquationSet_t", 1, NULL);
		cg_model_write("ThermalConductivityModel_t", CG_ConstantPrandtl);

		iterPrm = G_Domain.mapCasePrms_real.find("Pr");
		if (iterPrm != G_Domain.mapCasePrms_real.end())
		{
			const double Pr = iterPrm->second;
			cg_goto(f, iBase, "FlowEquationSet_t", 1, "ThermalConductivityModel_t", 1, NULL);
			cg_array_write("Prandtl", CG_RealDouble, 1, &one, &Pr);
			cg_goto(f, iBase, "FlowEquationSet_t", 1, "ThermalConductivityModel_t", 1, "DataArray_t", 1, NULL);
			cg_dataclass_write(CG_NondimensionalParameter);
		}
	}
	*/


	// Solution time
	{
		cg_biter_write(f, iBase, "TimeIterValues", 1/*time steps count*/);
		cg_goto(f, iBase, "BaseIterativeData_t", 1, NULL);

		const double t = G_State.time;

		const cgsize_t len = 1;
		cg_array_write("TimeValues", CG_RealDouble, 1, &len, &t);
		//cg_array_write("TimeStepValues"/*non-standart*/, CG_RealDouble, 1, &len, &dt);
	}

	return true;
}
