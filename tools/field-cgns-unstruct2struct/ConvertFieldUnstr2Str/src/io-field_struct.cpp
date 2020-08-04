///////////////////////////////////////////////////////////////////////////////
// Name:        io-field.cpp
// Purpose:     Input/Output functions
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////
#include "mpi.h"

#include "cgnslib_wrapper.h"
//#include "mem-usage.h"

#include "common_data_struct.h"

#include "common_procs_struct.h"

#include "logging.h"
#include "settings.h"
#include "io-field_struct.h"

//-----------------------------------------------------------------------------

static const char g_szCGBase[] = "HSFlow";
static const char* g_cgCoordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };
//-----------------------------------------------------------------------------

//
// Forward declarations
//
static double loadField_cgns(const std::string& fileName, const short time_layer);
static bool read_zone_cgns(const int fileID, const int iBase, const int idxZne,
	const int nxyz[], TpakArraysDyn<double>& grid, TpakArraysDyn<double>& field);
static double loadField_legacy( const std::string& fileName, const short time_layer);

static bool writeMetaInfoToCGNS(const int fileID, const int iBase, const short time_layer);
static std::string cgns_face_name(const int zone_idx, const TZoneFacePos face_pos);
//-----------------------------------------------------------------------------

//
// Helper classes
//

/**
 * Double or float array wrapper
 */
class TDFArray
{
	union {
		double* _D;
		float*  _F;
		void*   _raw;
	};
	bool is_double;

	TDFArray() = delete;
	TDFArray(TDFArray&) = delete;
	void operator=(TDFArray&) = delete;

public:
	TDFArray(bool d) : _raw(nullptr), is_double(d) { ; }

	void alloc(size_t size) {
		if( is_double )
			_D = new double[size];
		else
			_F = new float[size];
	}

	template<typename T>
	void set(size_t n, const T& val) {
		if( is_double )
			_D[n] = static_cast<double>(val);
		else
			_F[n] = static_cast<float>(val);
	}

	void* data() {
		return _raw;
	}

	/// Data type of stored data in CGNS & MPI terms
	struct TDataType {
		CG_DataType_t cgns;
		MPI_Datatype  mpi;
	};
	const TDataType type() const {
		if( is_double )
			return { CG_RealDouble, MPI_DOUBLE };
		else
			return { CG_RealSingle, MPI_FLOAT };
	}

	~TDFArray() {
		if( is_double )
			delete[] _D;
		else
			delete[] _F;
	}
};
//-----------------------------------------------------------------------------


/**
 *  Initializes field by freestream values or loading from file
 *
**/
bool initField()
{
	TLogSyncGuard logGuard;

	//
	// Memory allocation for field data
	//

	// Space to store field data of zones local to current MPI rank
	for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
	{
		TZone& zne = G_Domain.Zones[b];
		const int N = G_Domain.nu * zne.nodes_count();

		// Current (n+1) time layer
		zne.U = new double[N];

		// Previous (n) time layer
		if( ! zne.isFrozen )
			zne.Un = new double[N];

		// Pre-previous (n-1) time layer
		if( ! zne.isFrozen && G_TimePrm.approxOrder==2 )
			zne.Unm1 = new double[N];
	}

	hsLogDebug( "Field allocated. Memory: %.2f MiB", getCurrentRSS_MiB() );
	logGuard.Flush();

	double time = -1;

	//
	// Fill-in uniform field data
	//
	if( g_genOpts.strInitFieldFN.empty() )
	{
		const double* Uinf = G_Domain.phys->ref_state.UUinf;

		for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
		{
			TZone& zne = G_Domain.Zones[b];
			if( zne.isFrozen )
			{
				hsLogWarning(
					"Using uniform field for frozen zone '%s'#%d", zne.szName, b+1 );
			}

			for( int k = zne.ks; k <= zne.ke; ++k ){
			for( int j = zne.js; j <= zne.je; ++j ){
			for( int i = zne.is; i <= zne.ie; ++i )
			{
				const int pos = (zne.flatIdx(i,j,k) - 1)*G_Domain.nu;
				copy_n(Uinf, G_Domain.nu, zne.U + pos);
			}}}
		}

		time = 0.0;
	}

	//
	// Load initial estimation from file
	//
	else
	{
		hsLogMessage( "* Reading initial field '%s'...",
			g_genOpts.strInitFieldFN.c_str() );

		time = loadField(g_genOpts.strInitFieldFN, 0);
		if( time < 0 )
			return false;
	}

	// Set initial time in time-stepping
	G_Plugins.get_discret_caps().scheme->tab.tau =
		(g_genOpts.timeStart >= 0) ? g_genOpts.timeStart : time;

	return true;
}
//-----------------------------------------------------------------------------

/**
 * Loads field data from the file
 *
 * @param[in] fileName
 * @param[in] time_layer - time layer to load field into
 *                         0, 1 or 2 -> current, previous, pre-previous
 *
 * @retval >0  time value obtained from the file
 * @retval  0  if file has no time
 * @retval -1  failure
**/
double loadField(const std::string& fileName, const short time_layer)
{
	if( time_layer < 0 || time_layer > G_TimePrm.approxOrder ){
		hsLogError("Trying to load field into nonexistent time layer 'n%+d'", 1-time_layer);
		return -1;
	}

	double time = -1;
	if( hs_string_ends_with(fileName, ".cgns") )
		time = loadField_cgns(fileName, time_layer);
	else // .ttl, .hsx
		time = loadField_legacy(fileName, time_layer);

	MPI_Barrier( MPI_COMM_WORLD );
	hsLogWTime();
	hsLogMessage(" ");

	return time;
}
//-----------------------------------------------------------------------------


/**
 * Loads field data from the multiblock CGNS file
 *
 * @param[in] fileName
 * @param[in] time_layer - time layer to load field into
 *                         0, 1 or 2 -> current, previous, pre-previous
 *
 * @retval >0  time value obtained from the file
 * @retval  0  if file has no time
 * @retval -1  failure
**/
static double loadField_cgns( const std::string& fileName, const short time_layer)
{
	double time = -1.0;

	assert( time_layer >= 0 && time_layer <= G_TimePrm.approxOrder );
	TLogSyncGuard logGuard;

	const int &nu = G_Domain.nu;

	int f = -1;
	const int iBase = 1;  // assume only one base in the file

	char ok = 1;
if( G_State.mpiRank == 0 ) do
{
	ok = 0;

	if( cg_open( fileName.c_str(), CG_MODE_READ, &f ) != CG_OK )
	{
		hsLogError( "Can't open field file '%s' for reading: %s",
			fileName.c_str(), cg_get_error() );
		break;
	}

	char szName[33];

	// Space dimensions
	int dimCell = 0, dimPhys = 0;
	if( cg_base_read(f,iBase,  szName,&dimCell,&dimPhys) != CG_OK )
	{
		hsLogError( "Can't read CGNS base node from '%s' ( %s )",
			fileName.c_str(), cg_get_error() );
		break;
	}

	if( dimCell != G_Domain.nDim ) {
		hsLogError("CGNS: Inconsistent space dimensions");
		break;
	}

	// Number of zones (aka blocks)
	int nZones = 0;  cg_nzones(f,iBase, &nZones);
	if( nZones != G_Domain.nZones ) {
		hsLogError( "CGNS: Inconsistent number of zones" );
		break;
	}

	// Get time
	time = 0.0;
	do {
		int nTmSteps = 0;  cg_biter_read(f,iBase, szName, &nTmSteps);
		if( nTmSteps < 1 ) break;

		if( cg_goto(f,iBase,"BaseIterativeData_t",1, NULL) != CG_OK )
			break;

		int nArrs = 0;  cg_narrays(&nArrs);
		if( nArrs < 1 )  break;

		for( int ai = 1; ai <= nArrs; ++ai )
		{
			CG_DataType_t type;  int dim;  cgsize_t len[3];
			cg_array_info(ai, szName, &type, &dim, len);
			if( strcmp(szName, "TimeValues") == 0 ) {
				cg_array_read_as(ai, CG_RealDouble, &time);
				break;
			}
		}
	} while(false);

	ok = 1;
} while(false); // if( G_State.mpiRank == 0 )

	MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, PETSC_COMM_WORLD);
	if( ! ok )
		return -1;

	MPI_Bcast(&time, 1, MPI_DOUBLE, 0/*root*/, PETSC_COMM_WORLD);

	// Loop through zones
	for( int zi = 0; zi < G_Domain.nZones; ++zi )
	{
		TZone& zne = G_Domain.Zones[zi];
		if( zne.isFrozen && time_layer > 0 ) continue;
		double* dstU = zne.UU[time_layer];

		// Original zone size
		// (not accounting added ghosts & skipped layers)
		const int
			is0 = zne.is - (zne.Faces[faceXmin].isSkipped ? 1 : 0),
			ie0 = zne.ie + (zne.Faces[faceXmax].isSkipped ? 1 : 0),
			js0 = zne.js - (zne.Faces[faceYmin].isSkipped ? 1 : 0),
			je0 = zne.je + (zne.Faces[faceYmax].isSkipped ? 1 : 0),
			ks0 = zne.ks - (zne.Faces[faceZmin].isSkipped ? 1 : 0),
			ke0 = zne.ke + (zne.Faces[faceZmax].isSkipped ? 1 : 0);

		const int nxyz0[] = { ie0 - is0 + 1,  je0 - js0 + 1,  ke0 - ks0 + 1 };

		// Data of the zone excluding ghosts!!!
		TpakArraysDyn<double> newField, newGrid;
		if( (G_Domain.bs <= zi && zi <= G_Domain.be) || G_State.mpiRank == 0 )
		{
			//newGrid.reset(G_Domain.nDim, nxzy0[0]*nxzy0[1]*nxzy0[2]);
			newField.reset(G_Domain.nu, nxyz0[0]*nxyz0[1]*nxyz0[2]);
		}

		const int mpiTag = 'f'+'l'+'d' + zi;

		if( G_State.mpiRank == 0 )
		{
			if( ! read_zone_cgns(f,iBase,zi,nxyz0, newGrid, newField) )  ok = 0;

			// Send data from root to the zone's owner
			const int& rankDst = G_State.map_zone2rank[zi];
			if( rankDst != 0 )  // don't send to myself
				MPI_Ssend(newField.data(), (ok ? newField.size() : 0), MPI_DOUBLE, rankDst, mpiTag, PETSC_COMM_WORLD);
		}
		else if( G_Domain.bs <= zi && zi <= G_Domain.be )
		{
			MPI_Recv(newField.data(), newField.size(), MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Copy field to internal structs
		// NB: skipped grid layers are also filled to support immediate saving (e.g. slices) without scatternig ghosts
		if( G_Domain.bs <= zi && zi <= G_Domain.be )
		{
			for( int f = 0; f < nu; ++f ) // loop through field variables
			{
				const double* __restrict srcU = newField[f];
				for( int k = ks0;  k <= ke0;  ++k ){
				for( int j = js0;  j <= je0;  ++j ){
				for( int i = is0;  i <= ie0;  ++i )
				{
					dstU[nu*(zne.flatIdx(i,j,k) - 1) + f] = *srcU++;
				}}}
			}
		}
	} // loop through zones

	MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, PETSC_COMM_WORLD);
	if( ! ok )
		return -1;

	if( G_State.mpiRank == 0 )
		cg_close(f);

	return time;
}
//-----------------------------------------------------------------------------


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
static bool read_zone_cgns(const int fileID, const int iBase, const int idxZne,
	const int nxyz[], TpakArraysDyn<double>& grid, TpakArraysDyn<double>& field )
{
	const int& cgZneID = G_Domain.map_iZne2cgID[idxZne];

	// Get zone size and name
	char szZone[33];
	cgsize_t isize[9];  // NVertexI, NVertexJ, NVertexK,
	                    // NCellI, NCellJ, NCellK,
	                    // NBoundVertexI, NBoundVertexJ, NBoundVertexK
	if( cg_zone_read(fileID,iBase,cgZneID,  szZone,isize) != CG_OK ) {
		hsLogError( "Can't read CGNS zone #%d ( %s )",  cgZneID, cg_get_error() );
		return false;
	}

	if( strcmp(szZone, G_Domain.Zones[idxZne].szName) != 0 ) {
		hsLogWarning( "Inconsistent zone #%d names: '%s' -> '%s'",
			cgZneID,  szZone, G_Domain.Zones[idxZne].szName  );
	}

	// Block size without ghosts
	cgsize_t  nx = isize[0],  ny = isize[1],  nz = (G_Domain.nDim==3) ? isize[2] : 1;
	if( nxyz[0] != nx || nxyz[1] != ny || nxyz[2] != nz )
	{
		hsLogError( "Inconsistent zone '%s'#%d dimensions: %dx%dx%d -> %dx%dx%d",
			szZone, cgZneID,
			(int)nx, (int)ny, (int)nz,   nxyz[0], nxyz[1], nxyz[2]  );
		return false;
	}

	// Indexes faces
	cgsize_t irmin[3] = {1, 1, 1};
	cgsize_t irmax[3] = {nx, ny, nz};

	//
	// Read grid coordinates
	//
	if( grid.size() > 0 )
	for( int iCoord = 0; iCoord < G_Domain.nDim; ++iCoord )
	{
		const char* name = g_cgCoordNames[iCoord];
		if( cg_coord_read( fileID,iBase,cgZneID,name,CG_RealDouble, irmin,irmax,  grid[iCoord]) != CG_OK )
		{
			hsLogError( "Can't read %s from zone '%s'#%d ( %s )",
				name,   szZone,cgZneID,   cg_get_error()  );
			return false;
		}
	} // for iCoord


	//
	// Get solution info
	// FIXME: flow assumed existing
	//
	int iSol = 1;
	{
		CG_GridLocation_t loc;   char cgName[33];
		if( cg_sol_info(fileID,iBase,cgZneID,iSol,  cgName,&loc) != CG_OK )
		{
			hsLogError( "Can't read flow info from zone '%s'#%d ( %s )",
				szZone, cgZneID,   cg_get_error()  );
			return false;
		}

		if( loc != CG_Vertex )
		{
			hsLogError("CGNS: GridLocation must be Vertex");
			return false;
		}
	}

	//
	// Read functions
	//
	assert( field.size() > 0 );
	for( int iFun = 0; iFun < G_Domain.nu; ++iFun )
	{
		double* __restrict U = field[iFun];
		const char* fun_name = G_Domain.phys->vecFuncNames[iFun].c_str();

		int r = cg_field_read( fileID,iBase,cgZneID,iSol,fun_name,
			CG_RealDouble, irmin,irmax, U );

		if( r != CG_OK && r != CG_NODE_NOT_FOUND ) {
			hsLogError( "Can't read '%s' from zone '%s'#%d ( %s )",
				fun_name,   szZone, cgZneID,   cg_get_error()  );
			return false;
		}

		if( r == CG_NODE_NOT_FOUND ) {
			hsLogWarning( "%s not found in zone '%s'#%d, using free stream values",
				fun_name,   szZone, cgZneID  );

			std::fill_n(U, field.subsize(), G_Domain.phys->ref_state.UUinf[iFun] );
			continue;
		}

		/// Multi-line message builder
		class Tmlines {
			std::string text, line;

			void add_line() {
				text += "\n   " + line;
				line.clear();
			}

		public:
			Tmlines& operator<<(const std::string& msg) {
				line += msg + "; ";
				if( line.length() > 80 )  add_line();
				return *this;
			}

			const char* c_str() {
				if( ! line.empty() )  add_line();
				return text.c_str();
			}
		};

		// Filter negative values
		if( (g_genOpts.init_ltzero_filter > 0) &&
			( strcmp(fun_name, "Pressure")    == 0    ||
			  strcmp(fun_name, "Temperature") == 0    ||
			  strstr(fun_name, "TurbulentEnergy") == fun_name ) ) // begins with
		{
			const double fix = g_genOpts.init_ltzero_filter * G_Domain.phys->ref_state.UUinf[iFun];
			Tmlines msg;

			int m = 0;
			for( int k = 1; k <= nz; ++k ){
			for( int j = 1; j <= ny; ++j ){
			for( int i = 1; i <= nx; ++i )
			{
				double& val = U[m++];
				if( val < 0 ) {
					msg << hs_string_format("%.2g @(%d,%d,%d)", val,  i, j, k);
					val = fix;
				}
			}}}

			if( msg.c_str()[0] != '\0' )
				hsLogWarning( "Non-positive %s detected in '%s' zone: %s\n"
					"   Using %.3g (%.1f%% of freestream) instead...",
					fun_name, szZone,  msg.c_str(),
					fix, g_genOpts.init_ltzero_filter*100. );
		} // filer negative

		// Filter NaN
		{
			const double freestream_ratio = 0.5;
			const double fix = freestream_ratio * G_Domain.phys->ref_state.UUinf[iFun];
			Tmlines msg;

			int m = 0;
			for( int k = 1; k <= nz; ++k ){
			for( int j = 1; j <= ny; ++j ){
			for( int i = 1; i <= nx; ++i )
			{
				double& val = U[m++];
				if( isnan(val) || val > 1e66 || val < -1e66 ) {
					msg << hs_string_format("@(%d,%d,%d)",  i, j, k);
					val = fix;
				}
			}}}

			if( msg.c_str()[0] != '\0' )
				hsLogWarning( "NaN %s detected in '%s' zone: %s\n"
					"   Using %.3g (%.1f%% of freestream) instead...",
					fun_name, szZone,  msg.c_str(), fix, freestream_ratio*100. );
		} // filer NaN

	} //for iFun

	return true;
}
//-----------------------------------------------------------------------------


/**
 * Loads field data from the file in legacy format (.ttl, .hsx)
 *
 * @param[in] fileName
 * @param[in] time_layer - time layer to load field into
 *                         0, 1 or 2 -> current, previous, pre-previous
 *
 * @retval >0  time value obtained from the file
 * @retval  0  if file has no time
 * @retval -1  failure
**/
static double loadField_legacy(const std::string& fileName, const short time_layer)
{
	assert( time_layer >= 0 && time_layer <= G_TimePrm.approxOrder );
	if( G_Domain.nZones > 1 )
	{
		hsLogError( "Can't load multi-block field from legacy file format" );
		return -1;
	}

	const int iBlock = 0;
	TZone& blk = G_Domain.Zones[iBlock];
	double* dstU = blk.UU[time_layer];

	TDims dims;
	std::map<std::string, double> realParams;
	TdataArrayWrap funcArray(sizeof(double)),  gridArray(sizeof(double));

	if( ! gg_legacyFieldIO.Load( fileName,
			nullptr/*date*/, nullptr/*comments*/, nullptr/*func names*/,
			nullptr/*intPrms*/, &realParams,  dims, &funcArray, &gridArray) )
	{
		hsLogError( "%s", gg_legacyFieldIO.ErrMsg() );
		return -1;
	}

	if( dims.numFunc != G_Domain.nu )
	{
		hsLogError("Loaded field has incorrect number of functions");
		funcArray.freeMem();  gridArray.freeMem();
		return -1;
	}

	if( dims.numX != blk.nx || dims.numY != blk.ny || dims.numZ != blk.nz )
	{
		hsLogError(
			"Zone has different dimensions in the loaded field (%dx%dx%d) and in the current grid (%dx%dx%d).",
			dims.numX, dims.numY, dims.numZ,
			blk.nx, blk.ny, blk.nz
		);
		return -1;
	}

	memcpy(dstU, funcArray.A, G_Domain.nu * blk.nx * blk.ny * blk.nz * sizeof(double) );

	funcArray.freeMem();   gridArray.freeMem();

	return realParams["time"];
}
//-----------------------------------------------------------------------------


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
			   const short time_layer, bool isDouble)
{
	assert( time_layer >= 0 && time_layer <= G_TimePrm.approxOrder );

	int f = -1, fGrid = -1;  // file descriptors for field & grid
	std::string fn_grid;  // grid file name relative to the field
	int iBase = -1,  iBaseGrid = -1;  // CGNS base ids

	TLogSyncGuard logGuard;
	short ok = 1;  // error code for MPI broadcasting
if( G_State.mpiRank == 0 ) do
{
	ok = 0;

	const std::string path_field = g_CASE_RESULTS_DIR + fileName + ".cgns";

	if( cg_open( path_field.c_str(), CG_MODE_WRITE, &f ) != CG_OK )
	{
		hsLogError( "Can't open file '%s' for writing ( %s )",
			path_field.c_str(),   cg_get_error()  );
		break;
	}
	hsLogMessage( "\n* Saving field to file '%s'...", path_field.c_str() );

	fGrid = f;
	std::string path_grid = path_field;

	if( ! gridFileName.empty() )
	{
		//
		// Create separate grid file
		//
		fn_grid = gridFileName + ".cgns";
		path_grid = g_CASE_RESULTS_DIR + fn_grid;

		static bool needToSave = true;
		if( needToSave )
		{
			if( cg_open( path_grid.c_str(), CG_MODE_WRITE, &fGrid ) == CG_OK )
				needToSave = false;
			else {
				hsLogError( "Can't open file '%s' for writing ( %s )",
					path_grid.c_str(),   cg_get_error()  );
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
	if( cg_base_write(f,g_szCGBase, G_Domain.nDim,G_Domain.nDim, &iBase) != CG_OK )
	{
		hsLogError( "Can't write base into '%s' ( %s )",
			path_field.c_str(),   cg_get_error()  );
		break;
	}

	// Write Reference state, Equations info, time, etc
	if( ! writeMetaInfoToCGNS(f, iBase, time_layer) )
		break;

	iBaseGrid = iBase;
	if( fGrid != f && fGrid >= 0 ) // grid in separate file
	{
		if( cg_base_write(fGrid,g_szCGBase, G_Domain.nDim,G_Domain.nDim, &iBaseGrid) != CG_OK )
		{
			hsLogError( "Can't write CGNS base node into '%s' ( %s )",
				path_grid.c_str(),   cg_get_error()  );
			break;
		}
	}

	ok = 1;
} while(false); //if( G_State.mpiRank == 0 )

	MPI_Bcast(&ok, 1, MPI_SHORT, 0/*root*/, PETSC_COMM_WORLD);
	if( ! ok )
		return false;

	// Inform all ranks if we need grid coordinates or not
	MPI_Bcast(&fGrid, 1, MPI_INT, 0/*root*/, PETSC_COMM_WORLD);

	//
	// Sorted 1-based CGNS zone IDs and corresponding internal 0-based zone indices
	static int* map_cgID2iZne = nullptr;
	if( ! map_cgID2iZne )
	{
		map_cgID2iZne = new int[ G_Domain.nZones ];
		for( int b = 0; b < G_Domain.nZones; ++b )
			map_cgID2iZne[ G_Domain.map_iZne2cgID[b] - 1 ] = b;
	}

	//
	// Write field & grid data
	//
	for( int cgZneID = 1; cgZneID <= G_Domain.nZones; ++cgZneID )
	{
		const int& zi = map_cgID2iZne[cgZneID - 1];
		int iZone = -1;   int iZoneGrid = -1;

		TZone& zne = G_Domain.Zones[zi];

		// Indices including skipped grid layers
		const int is1 = zne.is - (zne.Faces[faceXmin].isSkipped ? 1 : 0);
		const int ie1 = zne.ie + (zne.Faces[faceXmax].isSkipped ? 1 : 0);
		const int js1 = zne.js - (zne.Faces[faceYmin].isSkipped ? 1 : 0);
		const int je1 = zne.je + (zne.Faces[faceYmax].isSkipped ? 1 : 0);
		const int ks1 = zne.ks - (zne.Faces[faceZmin].isSkipped ? 1 : 0);
		const int ke1 = zne.ke + (zne.Faces[faceZmax].isSkipped ? 1 : 0);

		// Size without ghost nodes
		const int nx0 = ie1 - is1 + 1;
		const int ny0 = je1 - js1 + 1;
		const int nz0 = ke1 - ks1 + 1;
		const int nxyz0 = nx0*ny0*nz0;

		if( G_State.mpiRank == 0 )
		{
			// Zone size packed in CGNS format
			cgsize_t isize[9];
			if( G_Domain.nDim == 2 )
			{
				isize[0] = nx0;    isize[1] = ny0;   // NVertexI, NVertexJ
				isize[2] = nx0-1;  isize[3] = ny0-1; // NCellI, NCellJ
				isize[4] = 0;      isize[5] = 0;     // NBoundVertexI, NBoundVertexJ
			}
			else
			{
				isize[0] = nx0;    isize[1] = ny0;    isize[2] = nz0;
				isize[3] = nx0-1;  isize[4] = ny0-1;  isize[5] = nz0-1;
				isize[6] = 0;      isize[7] = 0;      isize[8] = 0;
			}

			// Zone for fields
			cg_zone_write(f,iBase,zne.szName, isize,CG_Structured, &iZone);

			// Zone for grid coords
			iZoneGrid = iZone;
			if( fGrid != f && fGrid >= 0 ) // grid in separate file
			{
				cg_zone_write(fGrid,iBaseGrid,zne.szName,isize,CG_Structured, &iZoneGrid);
			}

			//
			// Make link to the grid in separate file
			//
			if( fGrid != f )
			{
				// CGNS-node path inside the separate grid file
				const char* label = "GridCoordinates";
				const std::string& path = hs_string_format("/%s/%s/%s", g_szCGBase, zne.szName, label);

				// Move CGNS file position to the current zone then write link
				cg_goto(f,iBase, "Zone_t",iZone, NULL);
				cg_link_write(label, fn_grid.c_str(), path.c_str() );
			}
		}


		TDFArray Vals(isDouble);
		if( (G_Domain.bs <= zi && zi <= G_Domain.be) || G_State.mpiRank == 0 )
			Vals.alloc(nxyz0);

		//
		// Grid coords
		//
		if( fGrid >= 0 )
		{
			for( int m = 0;  m < G_Domain.nDim; ++m )
			{
				const int mpiTag = 'g'+'r'+'d' + m;

				if( G_Domain.bs <= zi && zi <= G_Domain.be )
				{
					//
					// Worker mpi-rank -> fill-in and send
					//
					size_t c = 0;
					for( int k = ks1; k <= ke1; ++k ){
					for( int j = js1; j <= je1; ++j ){
					for( int i = is1; i <= ie1; ++i )
					{
						switch( 10*G_Domain.nDim + m )
						{
						// 2D
						case 20:  Vals.set(c, zne.coordIJ(i,j).x);     break;
						case 21:  Vals.set(c, zne.coordIJ(i,j).y);     break;

						// 3D
						case 30:  Vals.set(c, zne.coordIJK(i,j,k).x);  break;
						case 31:  Vals.set(c, zne.coordIJK(i,j,k).y);  break;
						case 32:  Vals.set(c, zne.coordIJK(i,j,k).z);  break;
						}
						++c;
					}}}

					if( G_State.mpiRank != 0 )
					{
						// NB: Don't use MPI_Send - it may flood the root MPI rank such that MPI_Recv fails
						MPI_Ssend(Vals.data(), nxyz0, Vals.type().mpi, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
					}
				}

				if( G_State.mpiRank == 0 )
				{
					//
					// Root mpi-rank -> receive and write
					//
					const int& rankSrc = G_State.map_zone2rank[zi];
					if( rankSrc != 0 ) // don't receive from myself
					{
						const int nn = ok ? nxyz0 : 0;  // if not OK, do a dummy recieve to unblock sender
						MPI_Recv( Vals.data(), nn, Vals.type().mpi,
							rankSrc, mpiTag, PETSC_COMM_WORLD,
							MPI_STATUS_IGNORE  // don't use NULL as status, MPI_Ssend may get stuck (i.e. in Intel MPI 5.0.1)
						);
					}

					int iCoord = -1;
					const char* name = g_cgCoordNames[m];
					if( ok )
					if( cg_coord_write(fGrid,iBaseGrid,iZoneGrid,  Vals.type().cgns,name, Vals.data(), &iCoord) != CG_OK )
					{
						hsLogError( "Can't write %s in zone %s#%d ( %s )",
							name,   zne.szName, iZoneGrid,   cg_get_error()  );
						ok = 0;  // don't return, continue recieving data from other ranks
					}
				}
			} // for( int m = 0;  m < G_Domain.nDim; ++m )
		} // if( fGrid >=0 )

		//
		// Wall BC info
		//
		for( int iFace = 0; iFace < 2*G_Domain.nDim; ++iFace )
		{
			const TZoneFace& face = zne.Faces[iFace];
			if( ! face.BC ) continue;
			if( ! face.BC->isWall() ) continue;

			const TZoneFacePos face_pos = static_cast<TZoneFacePos>(iFace);
			cgsize_t iRng[6] = { 0 };  // BC patch range

			int m = 0;
			switch( face_pos )
			{
			case faceXmin:
				iRng[m++] = 1;  iRng[m++] = 1;    if(G_Domain.nDim==3) iRng[m++] = 1;
				iRng[m++] = 1;  iRng[m++] = ny0;  if(G_Domain.nDim==3) iRng[m++] = nz0;
				break;

			case faceXmax:
				iRng[m++] = nx0;  iRng[m++] = 1;    if(G_Domain.nDim==3) iRng[m++] = 1;
				iRng[m++] = nx0;  iRng[m++] = ny0;  if(G_Domain.nDim==3) iRng[m++] = nz0;
				break;

			case faceYmin:
				iRng[m++] = 1;     iRng[m++] = 1;   if(G_Domain.nDim==3) iRng[m++] = 1;
				iRng[m++] = nx0;   iRng[m++] = 1;   if(G_Domain.nDim==3) iRng[m++] = nz0;
				break;

			case faceYmax:
				iRng[m++] = 1;     iRng[m++] = ny0;   if(G_Domain.nDim==3) iRng[m++] = 1;
				iRng[m++] = nx0;   iRng[m++] = ny0;   if(G_Domain.nDim==3) iRng[m++] = nz0;
				break;

			case faceZmin:
				iRng[m++] = 1;     iRng[m++] = 1;     if(G_Domain.nDim==3) iRng[m++] = 1;
				iRng[m++] = nx0;   iRng[m++] = ny0;   if(G_Domain.nDim==3) iRng[m++] = 1;
				break;

			case faceZmax:
				iRng[m++] = 1;     iRng[m++] = 1;     if(G_Domain.nDim==3) iRng[m++] = nz0;
				iRng[m++] = nx0;   iRng[m++] = ny0;   if(G_Domain.nDim==3) iRng[m++] = nz0;
				break;
			}

			int iBC = -1;
			const std::string& name = cgns_face_name(zi, face_pos);
			cg_boco_write(f,iBase,iZone, name.c_str(), CG_FamilySpecified, CG_PointRange, 2, iRng, &iBC);

			// Assign family name to the BC
			cg_goto(f,iBase, "Zone_t",iZone, "ZoneBC",0, "BC_t",iBC, NULL);
			cg_famname_write( face.szBC );

			// Create Family node
			int fam = 0;
			if( cg_family_write(f, iBase, face.szBC, &fam) == CG_OK ) // may already exist
			{
				// Make this Family the BC descriptor
				int famBC = 0;   cg_fambc_write(f, iBase, fam, "FamBC", CG_BCWall, &famBC);
			}
		} // loop through faces


		//
		// Field data
		//
		int iSol = -1;
		if( G_State.mpiRank == 0 )
			cg_sol_write(f,iBase,iZone,"FlowSolution",CG_Vertex, &iSol);

		const double* srcU = (! zne.isFrozen) ? zne.UU[time_layer] : zne.U;
		for( int fun = 0; fun < G_Domain.nu; ++fun )
		{
			const int mpiTag = 's'+'o'+'l' + fun;

			if( G_Domain.bs <= zi && zi <= G_Domain.be )
			{
				//
				// Worker mpi-rank -> fill-in and send
				//
				size_t c = 0;
				for( int k = ks1; k <= ke1; ++k ){
				for( int j = js1; j <= je1; ++j ){
				for( int i = is1; i <= ie1; ++i )
				{
					int pos = G_Domain.nu*( zne.flatIdx(i,j,k) - 1 ) + fun;
					Vals.set(c++, srcU[pos]);
				}}}

				if( G_State.mpiRank != 0 )
				{
					// NB: Don't use MPI_Send - it may flood the root MPI rank such that MPI_Recv fails
					MPI_Ssend(Vals.data(), nxyz0, Vals.type().mpi, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
				}
			}

			if( G_State.mpiRank == 0 )
			{
				//
				// Root mpi-rank -> receive and write
				//
				const int& rankSrc = G_State.map_zone2rank[zi];
				if( rankSrc != 0 ) // don't receive from myself
				{
					// if writing was failed previously, then do a dummy recieve to unblock sender
					MPI_Recv(Vals.data(), (ok ? nxyz0 : 0), Vals.type().mpi,
						rankSrc, mpiTag, PETSC_COMM_WORLD,
						MPI_STATUS_IGNORE  // don't use NULL as status, MPI_Ssend may get stuck (i.e. in Intel MPI 5.0.1)
					);
				}

				int iField = -1;
				const char* name = G_Domain.phys->vecFuncNames[fun].c_str();
				if( ok )
				if( cg_field_write(f,iBase,iZone,iSol,  Vals.type().cgns,name, Vals.data(), &iField ) != CG_OK )
				{
					hsLogError( "Can't write %s in zone %s#%d ( %s )",
						name,   zne.szName, iZone,   cg_get_error()  );
					ok = 0;  // don't return and continue recieving data from other ranks
				}
			}
		} // for( fun )

		// Vals[] is deleted here

	}  // Loop through zones

	if( G_State.mpiRank == 0 )
	{
		cg_close( f );

		if( fGrid != f && fGrid >= 0 )
			cg_close( fGrid );
	}

	MPI_Barrier( MPI_COMM_WORLD );
	hsLogWTime();

	return ok;
}
//-----------------------------------------------------------------------------


/**
 *   Wall faces info in the whole computation domain
 *
 *   NB: Function is UNUSED for now!!!
 */
void collectWallFaces()
{
	struct TWallFace
	{
		int iZne; // 0-based number of a domain zone, containing the face
		int pos;  // position of wall face in the zone
	};
	static TWallFace* wallFaces = NULL;

	if( wallFaces ) return;

	// Wall faces in the current MPI rank
	const int nWFacesLocMax = (G_Domain.be - G_Domain.bs + 1) * 6;
	TWallFace* wallFacesLoc = new TWallFace[nWFacesLocMax];
	int nWFacesLoc = 0;

	for( int iZne = G_Domain.bs; iZne <= G_Domain.be; ++iZne )
	{
		TZone& zne = G_Domain.Zones[iZne];
		for( int iFace = 0; iFace < 2*G_Domain.nDim; ++iFace )
		{
			if( ! zne.Faces[iFace].BC ) continue;
			if( ! zne.Faces[iFace].BC->isWall() ) continue;

			TWallFace& wf = wallFacesLoc[nWFacesLoc++];
			wf.iZne = iZne;
			wf.pos = iFace;
		}
	}

	// Distribute faces to all MPI ranks

	// How much to receive from each MPI rank
	int* recvcounts = new int[ G_State.mpiNProcs ];
	MPI_Allgather( &nWFacesLoc, 1,MPI_INT,  recvcounts, 1,MPI_INT,  PETSC_COMM_WORLD );

	// Total faces count
	int nWallFaces = 0;
	for( int r = 0; r < G_State.mpiNProcs; ++r )
		nWallFaces += recvcounts[r];

	// Exit if no wall faces at all
	if (nWallFaces == 0) {
		delete[] wallFacesLoc;
		return;
	}

	// Where put the data received from different MPI ranks
	wallFaces = new TWallFace[nWallFaces];

	// Displacements from the beginning of wallFaces[]
	int* displs = new int[ G_State.mpiNProcs ];
	displs[0] = 0;
	for( int r = 1; r < G_State.mpiNProcs; ++r )
		displs[r] = displs[r-1] + recvcounts[r-1];


	// Create custom mpi datatype for TWall::TFace
	MPI_Datatype MPI_TwallFace;
	{
		int          counts[2]  = {    1,       1,   };
		MPI_Datatype types[2]   = { MPI_INT, MPI_INT };
		MPI_Aint     offsets[2] = { offsetof(TWallFace, iZne),
									offsetof(TWallFace, pos)   };
		MPI_Type_create_struct(2, counts, offsets, types, &MPI_TwallFace);
	}
	MPI_Type_commit(&MPI_TwallFace);


	// Do send/receive wall faces info
	MPI_Allgatherv(
		wallFacesLoc, nWFacesLoc, MPI_TwallFace,
		wallFaces, recvcounts, displs, MPI_TwallFace,
		PETSC_COMM_WORLD
	);

	delete[] wallFacesLoc;
	delete[] recvcounts;
	delete[] displs;
	MPI_Type_free(&MPI_TwallFace);
}
//-----------------------------------------------------------------------------


/**
 *  Saves field data on a wall to file
 *
 * @param[in] fileName - file name to save wall data (without .cgns extension)
 * @param[in] gridFileName - file name for grid, if NULL embed grid into field file
 *
**/
bool saveFieldWall(const std::string& fileName, const std::string& gridFileName)
{
	// directory to save fields to
	const std::string dir = std::string() + g_CASE_RESULTS_DIR + "wall/";

	if( G_Domain.nDim < 3 )
		return false;

	// SAVING ONLY WALLS (2D cells with 3D coordinates)
	const int cell_dim = 2;
	const int phys_dim = G_Domain.nDim;

	TLogSyncGuard logGuard;
	short ok = 1;  // assumed to be '0' or '1' exactly

	//
	// Open CGNS file
	//
	int f = -1;  int fGrid = -1;
	int iBase = -1;  int iBaseGrid = -1;
	std::string fnGrid;

	if( G_State.mpiRank == 0 ) do
	{
		ok = 0;

		// Path to field file
		const std::string fnField = dir + fileName + ".cgns";
		hs_dir_create(dir);  // NB: g_CASE_RESULTS_DIR should be already created

		const char* field_fn = fnField.c_str();

		if( cg_open( field_fn, CG_MODE_WRITE, &f ) != CG_OK )
		{
			hsLogError( "Can't open file '%s' for writing ( %s )", field_fn, cg_get_error() );
			ok = 0;  break;
		}

		hsLogMessage( "\n* Saving wall to '%s'...", field_fn );

		fGrid = f;
		fnGrid = fnField;
		if( ! gridFileName.empty() )  // Create separate grid file
		{
			// File name of grid
			fnGrid = dir + gridFileName + ".cgns";
			const char* grid_fn = fnGrid.c_str();

			static bool isGridNeeded = true;
			if( isGridNeeded )
			{
				if( cg_open( grid_fn, CG_MODE_WRITE, &fGrid ) != CG_OK )
				{
					hsLogError( "Can't open file '%s' for writing ( %s )",
						grid_fn, cg_get_error() );
					ok = 0;  break;
				}
				isGridNeeded = false;
			}
			else
				fGrid = -1;

			// Grid file name relative to field file
			fnGrid = gridFileName + ".cgns";
		}

		// Create base
		cg_base_write(f,g_szCGBase, cell_dim, phys_dim, &iBase);

		// Write Reference state, Equations info, time, etc
		writeMetaInfoToCGNS(f, iBase, 0);

		iBaseGrid = iBase;
		if( fGrid != f && fGrid >= 0 ) // grid in separate file
			cg_base_write(fGrid,g_szCGBase, cell_dim,phys_dim, &iBaseGrid);

		ok = 1;
	} while(false); //if( G_State.mpiRank == 0 )

	MPI_Bcast(&ok, 1, MPI_SHORT, 0/*root*/, PETSC_COMM_WORLD);
	if( ! ok )
		return false;

	// Inform all ranks if we need grid coordinates or not
	MPI_Bcast(&fGrid, 1, MPI_INT, 0/*root*/, PETSC_COMM_WORLD);

//---

	//
	// Loop through all zones looking for walls
	//
	for( int iZne = 0; iZne < G_Domain.nZones; ++iZne )
	{
		TZone& zne = G_Domain.Zones[iZne];
		for( int iFace = 0; iFace < 2*G_Domain.nDim; ++iFace )
		{
			const TZoneFace& face = zne.Faces[iFace];
			if( ! face.BC ) continue;
			if( ! face.BC->isWall() ) continue;

			const TZoneFacePos posWFace = static_cast<TZoneFacePos>(iFace);

			// Indices of grid nodes including skipped grid layers
			int is1 = zne.is - (zne.Faces[faceXmin].isSkipped ? 1 : 0);
			int ie1 = zne.ie + (zne.Faces[faceXmax].isSkipped ? 1 : 0);
			int js1 = zne.js - (zne.Faces[faceYmin].isSkipped ? 1 : 0);
			int je1 = zne.je + (zne.Faces[faceYmax].isSkipped ? 1 : 0);
			int ks1 = zne.ks - (zne.Faces[faceZmin].isSkipped ? 1 : 0);
			int ke1 = zne.ke + (zne.Faces[faceZmax].isSkipped ? 1 : 0);

			// Dimensions without ghost nodes
			const int nx0 = ie1 - is1 + 1;
			const int ny0 = je1 - js1 + 1;
			const int nz0 = ke1 - ks1 + 1;

			// Wall zone dimensions
			int n1, n2;
			switch( posWFace )
			{
			case faceXmin:   ie1 = is1;  n1 = ny0;  n2 = nz0;  break;
			case faceXmax:   is1 = ie1;  n1 = ny0;  n2 = nz0;  break;

			case faceYmin:   je1 = js1;  n1 = nx0;  n2 = nz0;  break;
			case faceYmax:   js1 = je1;  n1 = nx0;  n2 = nz0;  break;

			case faceZmin:   ke1 = ks1;  n1 = nx0;  n2 = ny0;  break;
			case faceZmax:   ks1 = ke1;  n1 = nx0;  n2 = ny0;  break;
			}


			int cgWZne = -1, cgWZneGrid = -1;  // CGNS IDs of wall zone & its grid
			if( G_State.mpiRank == 0 )
			{
				const std::string& cgName = cgns_face_name(iZne, posWFace);

				// Zone size packed in CGNS format
				cgsize_t isize[6];
				isize[0] = n1;      isize[1] = n2;     // NVertexI, NVertexJ
				isize[2] = n1 - 1;  isize[3] = n2 - 1; // NCellI, NCellJ
				isize[4] = 0;       isize[5] = 0;      // NBoundVertexI, NBoundVertexJ

				// CGNS zone for wall face
				cg_zone_write(f, iBase, cgName.c_str(), isize, CG_Structured, &cgWZne);

				// Zone for grid coords
				cgWZneGrid = cgWZne;
				if( fGrid != f && fGrid >= 0 ) // grid in separate file
					cg_zone_write(fGrid, iBaseGrid, cgName.c_str(), isize, CG_Structured, &cgWZneGrid);

				//
				// Make link to the grid in separate file
				//
				if( fGrid != f )
				{
					// CGNS-node path inside the separate grid file
					const char* label = "GridCoordinates";
					const std::string& path = hs_string_format("/%s/%s/%s", g_szCGBase, cgName.c_str(), label);

					// Move CGNS file position to the current zone then write link
					cg_goto(f, iBase, "Zone_t", cgWZne, NULL);
					cg_link_write(label, fnGrid.c_str(), path.c_str());
				}
			}

			double* Vals = new double[n1 * n2];

			//
			// Grid coords
			//
			if( fGrid >= 0 )
			{
				for( int m = 0; m < phys_dim; ++m )
				{
					const int mpiTag = 'w' + 'g' + 'r' + 'd' + m;

					if( G_Domain.bs <= iZne && iZne <= G_Domain.be )
					{
						//
						// Worker MPI rank -> fill-in and send
						//

						int cnt = 0;
						for( int k = ks1; k <= ke1; ++k ){
						for( int j = js1; j <= je1; ++j ){
						for( int i = is1; i <= ie1; ++i )
						{
							double& val = Vals[cnt++];
							switch( 10 * G_Domain.nDim + m )
							{
							//case 20:  val = zne.coordIJ(i,j).x;  break;  // 2D, coordX
							//case 21:  val = zne.coordIJ(i,j).y;  break;  // 2D, coordY
							case 30:  val = zne.coordIJK(i,j,k).x;  break;  // 3D, coordX
							case 31:  val = zne.coordIJK(i,j,k).y;  break;  // 3D, coordY
							case 32:  val = zne.coordIJK(i,j,k).z;  break;  // 3D, coordZ
							}
						}}}

						if( G_State.mpiRank != 0 )
							MPI_Ssend(Vals, n1*n2, MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
					}

					if( G_State.mpiRank == 0 )
					{
						//
						// Root MPI rank -> receive and write
						//
						const int& rankSrc = G_State.map_zone2rank[iZne];
						if( rankSrc != 0 ) // don't receive from myself
						{
							MPI_Recv(Vals, ok* n1*n2, MPI_DOUBLE, // dummy recieve if ok==0 to unblock sender
								rankSrc, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
						}

						int iCoord = -1;
						const char* name = g_cgCoordNames[m];
						if( ok )
						if( cg_coord_write(fGrid, iBaseGrid, cgWZneGrid, CG_RealDouble, name, Vals, &iCoord) != CG_OK )
						{
							hsLogError("Can't write %s to zone #%d ( %s )", name, cgWZneGrid, cg_get_error());
							ok = 0;  // don't return, continue recieving data from other ranks
						}
					}
				} // for( int m = 0;  m < G_Domain.nDim; ++m )
			} // if( fGrid >=0 )


			//
			// Field data
			//
			int iSol = -1;
			if( G_State.mpiRank == 0 )
				cg_sol_write(f, iBase, cgWZne, "FlowSolution", CG_Vertex, &iSol);

			const std::set<std::string>& setFuncs = g_genOpts.setWallFuncs;
			unsigned short n_saved_wall_funcs = 0;

			for( int fun = 0; fun < G_Domain.nu; ++fun )
			{
				const char* szFunName = G_Domain.phys->vecFuncNames[fun].c_str();
				if( setFuncs.count(szFunName) == 0 ) continue;
				++n_saved_wall_funcs;

				const int mpiTag = 'w'+'s'+'o'+'l' + fun;

				if( G_Domain.bs <= iZne && iZne <= G_Domain.be )
				{
					//
					// Worker MPI rank -> fill-in and send
					//
					int cnt = 0;
					for( int k = ks1; k <= ke1; ++k ){
					for( int j = js1; j <= je1; ++j ){
					for( int i = is1; i <= ie1; ++i )
					{
						const int pos = G_Domain.nu*(zne.flatIdx(i, j, k) - 1) + fun;
						Vals[cnt++] = zne.U[pos];
					}}}

					if( G_State.mpiRank != 0 )
						MPI_Ssend(Vals, n1*n2, MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
				}

				if( G_State.mpiRank == 0 )
				{
					//
					// Root MPI rank -> receive and write
					//
					const int& rankSrc = G_State.map_zone2rank[iZne];
					if( rankSrc != 0 )
						MPI_Recv(Vals, ok* n1*n2, MPI_DOUBLE,  // dummy recieve if ok==0 to unblock sender
							rankSrc, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

					int iField = -1;
					if( ok )
					if( cg_field_write(f,iBase,cgWZne, iSol, CG_RealDouble, szFunName,Vals, &iField) != CG_OK )
					{
						hsLogError( "Can't write %s to zone #%d ( %s )",
							szFunName, cgWZne,  cg_get_error() );
						ok = 0;  // don't return, continue recieving data from other ranks
					}
				}
			} // for( fun )


			//
			// Saving derived values on the wall
			//
			{
				enum { dv_CfX=0, dv_dTdn, dv_mu, dv_MAX };
				static const char* namesDrvdVals[] = { "CoefSkinFrictionX", "TemperatureGradientNormal", "ViscosityMolecular" };

				// Function index in solution array
				static int idxVx = -1, idxT = -1;
				if( idxVx < 0 )
				{
					for( int m = 0; m < G_Domain.nu; ++m )
					{
						if( G_Domain.phys->vecFuncNames[m] == "VelocityX" )
							idxVx = m;
						else if( G_Domain.phys->vecFuncNames[m] == "Temperature" )
							idxT = m;
					}
				}

				// Shift along grid line in wall normal direction
				int di = 0, dj = 0, dk = 0;
				switch( posWFace ) 
				{
				case faceXmin:   di = +1;   break;
				case faceXmax:   di = -1;   break;
				case faceYmin:   dj = +1;   break;
				case faceYmax:   dj = -1;   break;
				case faceZmin:   dk = +1;   break;
				case faceZmax:   dk = -1;   break;
				}

				for( int fun = 0; fun < dv_MAX; ++fun )
				{
					if( setFuncs.count(namesDrvdVals[fun]) == 0 )  continue;
					++n_saved_wall_funcs;

					const int mpiTag = 'w' + 's' + 'o' + 'l' + 'D' + fun;

					if( G_Domain.bs <= iZne && iZne <= G_Domain.be )
					{
						//
						// Worker MPI rank -> fill-in and send
						//
						int cnt = 0;
						for( int k = ks1; k <= ke1; ++k ){
						for( int j = js1; j <= je1; ++j ){
						for( int i = is1; i <= ie1; ++i )
						{
							// 1st, 2nd & 3rd nodes in normal-to-wall direction
							const double* const F1 = zne.U + (zne.flatIdx(i,     j,     k)      - 1)*G_Domain.nu;
							const double* const F2 = zne.U + (zne.flatIdx(i+di,  j+dj,  k+dk)   - 1)*G_Domain.nu;
							const double* const F3 = zne.U + (zne.flatIdx(i+2*di,j+2*dj,k+2*dk) - 1)*G_Domain.nu;

							const Tcoord3D
								&p1 = zne.coordIJK(i,   j,   k),
								&p2 = zne.coordIJK(i+di,j+dj,k+dk);
							const double h = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));

							switch( fun )
							{
							case dv_CfX:
								Vals[cnt++] = G_Domain.pfunViscosity(F1[idxT])/G_Domain.phys->ref_state.Re * (4*F2[idxVx] - 3*F1[idxVx] - F3[idxVx])/h;
								break;

							case dv_dTdn:
								Vals[cnt++] = (4*F2[idxT] - 3*F1[idxT] - F3[idxT])/h;
								break;

							case dv_mu:
								Vals[cnt++] = G_Domain.pfunViscosity(F1[idxT]);
								break;
							}
						}}}

						if( G_State.mpiRank != 0 )
							MPI_Ssend(Vals, n1*n2, MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
					}

					if( G_State.mpiRank == 0 )
					{
						//
						// Root MPI rank -> receive and write
						//
						const int& rankSrc = G_State.map_zone2rank[iZne];
						if( rankSrc != 0 )
							MPI_Recv(Vals, ok* n1*n2, MPI_DOUBLE, rankSrc, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

						int iField = -1;

						const char* name = namesDrvdVals[fun];
						if( ok )
						if( cg_field_write(f, iBase, cgWZne, iSol, CG_RealDouble, name, Vals, &iField) != CG_OK )
						{
							hsLogError("Can't write %s to zone #%d ( %s )", name, cgWZne, cg_get_error());
							ok = 0;
						}
					}
				}  // for (fun)
			} // derived values

			if( n_saved_wall_funcs < setFuncs.size() )
			{
				static bool err_reported = false;
				if( ! err_reported ){
					hsLogWarning("Some functions requested to save on the wall are not available");
					err_reported = true;
				}
			}

			delete[] Vals;
		} // loop through faces
	} // loop through zones

	if( G_State.mpiRank == 0 )
	{
		cg_close( f );

		if( fGrid != f && fGrid >= 0 )
			cg_close( fGrid );
	}

	MPI_Barrier( MPI_COMM_WORLD );
	hsLogWTime();

	return ok;
}
//-----------------------------------------------------------------------------


/**
 *  Saves field data at grid slices
 *  Info on the slices is stored in g_Slices
**/
bool saveFieldSlices(const std::string& baseFileName, const std::string& gridFileName)
{
	static bool* needToSaveGrd = NULL;
	if( ! needToSaveGrd )
	{
		needToSaveGrd = new bool[ g_Slices.size() ];
		for( int m = 0; m < g_Slices.size(); ++m )
			needToSaveGrd[m] = true;
	}

	hsLogMessage( "\n* Saving slices to files:" );
	short ok = 1;  // assumed to be '0' or '1' exactly
	TLogSyncGuard logGuard;

	// NB: g_Slices have been filled in doLoadGrid_cgns(...)
	std::vector<TSlice>::const_iterator itSlice = g_Slices.begin();
	int numSlice = 0;
	for( ; itSlice != g_Slices.end(); ++itSlice, ++numSlice )
	{
		const TSlice& slice = *itSlice;

		// Grid slice: 2D cells with 3D coordinates
		const int cell_dim = 2;
		const int phys_dim = G_Domain.nDim;

		int f = -1;  int fGrid = -1;
		int iBase = -1;  int iBaseGrid = -1;

		const std::string dir = g_CASE_RESULTS_DIR + slice.name;  // directory to save slices to
		std::string fnGrid;

		if( G_State.mpiRank == 0 ) do
		{
			// Path to field file
			std::string fnField = dir + "/" + baseFileName + ".cgns";
			hs_dir_create(dir);   // NB: strCASE_RESULTS_DIR should be already created

			const char* field_fn = fnField.c_str();

			if( cg_open( field_fn,CG_MODE_WRITE, &f ) != CG_OK )
			{
				hsLogError( "Can't open file '%s' for writing ( %s )",
					field_fn,   cg_get_error()   );
				ok = 0;  break;
			}
			hsLogMessage( "    %s", field_fn );

			fGrid = f;
			fnGrid = fnField;
			if( ! gridFileName.empty() )
			{
				// Create separate grid file

				// File name of grid
				fnGrid = dir + "/" + gridFileName + ".cgns";
				const char* grid_fn = fnGrid.c_str();

				if( needToSaveGrd[numSlice] )
				{
					if( cg_open( grid_fn, CG_MODE_WRITE, &fGrid ) != CG_OK )
					{
						hsLogError( "Can't open file '%s' for writing ( %s )",
							grid_fn, cg_get_error() );
						ok = 0;  break;
					}
					needToSaveGrd[numSlice] = false;
				}
				else
					fGrid = -1;

				// Grid file name relative to field file
				fnGrid = gridFileName + ".cgns";
			}

			// Create base
			cg_base_write(f,g_szCGBase, cell_dim, phys_dim, &iBase);

			// Write Reference state, Equations info, time, etc
			writeMetaInfoToCGNS(f, iBase, 0);

			iBaseGrid = iBase;
			if( fGrid != f && fGrid >= 0 ) // grid in separate file
			{
				cg_base_write(fGrid,g_szCGBase, cell_dim,phys_dim, &iBaseGrid);
			}
		} while(false); //if( G_State.mpiRank == 0 )

		MPI_Bcast(&ok, 1, MPI_SHORT, 0/*root*/, PETSC_COMM_WORLD);
		if( ! ok )
			return false;

		//---
		MPI_Bcast(&fGrid, 1, MPI_INT, 0/*root*/, PETSC_COMM_WORLD);

		//
		//  Looping through zones in current slice
		//
		std::vector<TSliceZneIdx>::const_iterator itSliceZne = slice.zn.begin();
		for( ; itSliceZne < slice.zn.end(); ++itSliceZne )
		{
 			const TSliceZneIdx& sliceZne = *itSliceZne;

			const int& b = sliceZne.iZne; // number of current zone
 			TZone& zne = G_Domain.Zones[b];

			const int
				&is = sliceZne.is, &ie = sliceZne.ie,
				&js = sliceZne.js, &je = sliceZne.je,
				&ks = sliceZne.ks, &ke = sliceZne.ke;

			int n1 = 0, n2 = 0;
			{
				const int nx = ie - is + 1;
				const int ny = je - js + 1;
				const int nz = ke - ks + 1;

				     if( nx == 1 ){   n1 = ny;   n2 = nz;   }
				else if( ny == 1 ){   n1 = nx;   n2 = nz;   }
				else if( nz == 1 ){   n1 = nx;   n2 = ny;   }
			}

			int iZone = -1;
			int iZoneGrid = -1;

			if( G_State.mpiRank == 0 )
			{
				// Zone size packed in CGNS format
				cgsize_t isize[6];
				isize[0] = n1;    isize[1] = n2;   // NVertexI, NVertexJ
				isize[2] = n1-1;  isize[3] = n2-1; // NCellI, NCellJ
				isize[4] = 0;     isize[5] = 0;    // NBoundVertexI, NBoundVertexJ

				// Zone for field data
				cg_zone_write(f,iBase,zne.szName, isize,CG_Structured, &iZone);

				// Zone for grid coords
				iZoneGrid = iZone;
				if( fGrid != f && fGrid >= 0 ) // grid in separate file
				{
					cg_zone_write(fGrid,iBaseGrid,zne.szName,isize,CG_Structured, &iZoneGrid);
				}

				//
				// Make link to the grid in separate file
				//
				if( fGrid != f )
				{
					// CGNS-node path inside the separate grid file
					const char* label = "GridCoordinates";
					const std::string& path = hs_string_format("/%s/%s/%s", g_szCGBase, zne.szName, label);

					// Move CGNS file position to the current zone then write link
					cg_goto(f,iBase, "Zone_t",iZone, NULL);
					cg_link_write(label, fnGrid.c_str(), path.c_str() );
				}
			}
			double* Vals = new double[ n1*n2 ];

			//
			// Grid coords
			//
			if( fGrid >= 0 ){
			for( int m = 0;  m < phys_dim; ++m )
			{
				const int mpiTag = 'g'+'r'+'d' + m;

				if( G_Domain.bs <= b && b <= G_Domain.be )
				{
					//
					// Worker #mpi-rank -> fill-in and send
					//
					int cnt = 0;
					for( int k = ks; k <= ke; ++k ){
					for( int j = js; j <= je; ++j ){
					for( int i = is; i <= ie; ++i )
					{
						double& val = Vals[ cnt++ ];
						switch( 10*G_Domain.nDim + m )
						{
						//case 20:  val = zne.coordIJ(i,j).x;  break;  // 2D, coordX
						//case 21:  val = zne.coordIJ(i,j).y;  break;  // 2D, coordY
						case 30:  val = zne.coordIJK(i,j,k).x;  break;  // 3D, coordX
						case 31:  val = zne.coordIJK(i,j,k).y;  break;  // 3D, coordY
						case 32:  val = zne.coordIJK(i,j,k).z;  break;  // 3D, coordZ
						}
					}}}

					if( G_State.mpiRank != 0 )
						// NB: Don't use MPI_Send - it may flood the root MPI rank such that MPI_Recv fails
						MPI_Ssend(Vals, n1*n2, MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
				}

				if( G_State.mpiRank == 0 )
				{
					//
					// Root mpi-rank -> receive and write
					//
					const int& rankSrc = G_State.map_zone2rank[b];
					if( rankSrc != 0 )
						MPI_Recv(Vals, ok* n1*n2, MPI_DOUBLE,  // dummy recieve if ok==0
							rankSrc, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

					int iCoord = -1;
					const char* name = g_cgCoordNames[m];
					if( cg_coord_write(fGrid,iBaseGrid,iZoneGrid, CG_RealDouble,name, Vals, &iCoord) != CG_OK )
					{
						hsLogError( "Can't write %s in zone %s#%d ( %s )",
							name,   zne.szName, iZoneGrid,   cg_get_error()  );
						ok = 0;  // don't return, just continue recieving from other ranks
					}
				}
			} // for( int m = 0;  m < G_Domain.nDim; ++m )
			} // if( fGrid >=0 )

			//
			// Field data
			//
			int iSol = -1;
			if( G_State.mpiRank == 0 )
				cg_sol_write(f,iBase,iZone,"FlowSolution",CG_Vertex, &iSol);

			///////////////////////////////////////
			// Saving flow data at current slice //
			///////////////////////////////////////
			for( int fun = 0; fun < G_Domain.nu; ++fun )
			{
				const int mpiTag = 's'+'o'+'l' + fun;

				if( G_Domain.bs <= b && b <= G_Domain.be )
				{
					//
					// Worker mpi-rank -> fill-in and send
					//
					int cnt = 0;
					for( int k = ks; k <= ke; ++k ){
					for( int j = js; j <= je; ++j ){
					for( int i = is; i <= ie; ++i )
					{
						int pos = G_Domain.nu*( zne.flatIdx(i,j,k) - 1 ) + fun;
						Vals[ cnt++ ] = zne.U[ pos ];
					}}}

					if( G_State.mpiRank != 0 )
						// NB: Don't use MPI_Send - it may flood the root MPI rank such that MPI_Recv fails
						MPI_Ssend(Vals, n1*n2, MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD);
				}

				if( G_State.mpiRank == 0 )
				{
					//
					// Root mpi-rank -> receive and write
					//
					const int& rankSrc = G_State.map_zone2rank[b];
					if( rankSrc != 0 )
						MPI_Recv(Vals, ok* n1*n2, MPI_DOUBLE,  // dummy recieve if ok==0
							rankSrc, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

					int iField = -1;
					const char* name = G_Domain.phys->vecFuncNames[fun].c_str();
					if( cg_field_write(f,iBase,iZone,iSol,   CG_RealDouble,name, Vals, &iField ) != CG_OK )
					{
						hsLogError(  "Can't write %s in zone %s#%d ( %s )",
							name,   zne.szName, iZone,   cg_get_error()   );
						ok = 0;  // don't return, continue recieving from other ranks
					}
				}
			} // for( fun )

			delete[] Vals;
		}  // Loop through zones in the slice

		if( G_State.mpiRank == 0 )
		{
			cg_close( f );

			if( fGrid != f && fGrid >= 0 )
				cg_close( fGrid );
		}
	} // loop through slices

	MPI_Barrier( MPI_COMM_WORLD );
	hsLogWTime();

	return ok;
}

//-----------------------------------------------------------------------------

/**
 *  Saves field data to file. Supports only single-block
 *
 * @param[in] aFileName - file name to save field data, if NULL generate from time
**/
bool saveField_legacy(const char* aFileName)
{
	// Saving from root rank only
	if( G_State.mpiRank != 0 )  return true;

	if( G_Domain.nZones > 1 )
	{
		hsLogWarning("Saving multiblock field using legacy file format. Saving the 1st block only...");
	}
	TZone& blk = G_Domain.Zones[0];
	const int &nx = blk.nx,  &ny = blk.ny,  &nz = blk.nz;

	const int &nu = G_Domain.nu,  &ndim = G_Domain.nDim;

	static double* funcArray = nullptr;
	static TdataArrayWrap funcWrapArray(nu * nx*ny*nz, funcArray);

	double* gridArray = nullptr;
	static TdataArrayWrap gridWrapArray(ndim * nx*ny*nz, gridArray);

	static std::map<std::string, int>  mapCasePrms_int;

	// Инициализация параметров для записи
	static bool isInit = true;
	if(isInit)
	{
		//Резервирование памяти для массивов данных в нужном для записи формате
		funcArray = (double*)funcWrapArray.allocMem();
		gridArray = (double*)gridWrapArray.allocMem();

		mapCasePrms_int["timeOrd"] = G_TimePrm.approxOrder;

		G_Domain.mapCasePrms_real["M"] = G_Domain.phys->ref_state.M;
		G_Domain.mapCasePrms_real["Re"] = G_Domain.phys->ref_state.Re;
		G_Domain.mapCasePrms_real["Tinf[K]"] = G_Domain.phys->ref_state.Tinf;
	}


	//Преобразование данных из рабочего формата в формат пригодный для сохранения:
	int laXY = 0,  laF = 0;
	for(int k=1; k<=nz; k++)
	for(int j=1; j<=ny; j++) //строки
	for(int i=1; i<=nx; i++) //столбцы
	{
		//Узел (#x,#y)=(i,j):
		if(gridArray)
		{
			if(nz<=1)
			{
				const Tcoord2D& p = blk.coordIJ(i,j);
				gridArray[laXY++] = p.x;
				gridArray[laXY++] = p.y;
			}
			else
			{
				const Tcoord3D& p = blk.coordIJK(i,j,k);
				gridArray[laXY++] = p.x;
				gridArray[laXY++] = p.y;
				gridArray[laXY++] = p.z;
			}
		} //физические X и Y

		for( int m=0;  m<nu;  ++m )
		{
			funcArray[laF] = blk.U[laF];
			++laF;
		}
	}

	// Получение полного времени:
	double time = G_Plugins.get_discret_caps().scheme->tab.tau;
	G_Domain.mapCasePrms_real["time"] = time;

	std::string fileName = g_CASE_RESULTS_DIR;
	if(aFileName)
		fileName += aFileName;
	else
		fileName += hs_string_format("%08.5f", time);

	std::string	strFunctionNames;
	for( int f = 0; f < nu; ++f )
		strFunctionNames += G_Domain.phys->vecFuncNames[f] + "\n";

	bool res = gg_legacyFieldIO.Save(
		fileName, ""/*comments*/, strFunctionNames, wfGZBin,
	    &mapCasePrms_int, &G_Domain.mapCasePrms_real,  TDims(nu, nx, ny, nz),
		&funcWrapArray,  &gridWrapArray,  "_grid"
	);
	if( ! res )
	{
		hsLogError( "%s", gg_legacyFieldIO.ErrMsg() );
		return false;
	}

	if(isInit){  gridWrapArray.freeMem();  isInit = false;  }

	hsLogMessage(  "* Field saved to file '%s'",  fileName.c_str()  );
	return true;
}
//-----------------------------------------------------------------------------


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
	cg_goto(f,iBase, NULL);
	cg_dataclass_write( CG_NormalizedByUnknownDimensional );

	//
	// Reference state
	//
	cg_goto(f,iBase, NULL);
	cg_state_write("ReferenceQuantities");
	{
		const cgsize_t one = 1;
		int count = 0;  // reference value number

		const hsflow::TFreestream& refQ = G_Domain.phys->ref_state;

		// Mach
		if( ! isnan(refQ.M) )
		{
			cg_goto(f,iBase,"ReferenceState_t",1, NULL);
			cg_array_write("Mach", CG_RealDouble,1,&one, &refQ.M);
			cg_goto(f,iBase,"ReferenceState_t",1, "DataArray_t",++count,NULL);
			cg_dataclass_write( CG_NondimensionalParameter );
		}

		// Re
		if( ! isnan(refQ.Re) )
		{
			cg_goto(f,iBase,"ReferenceState_t",1, NULL);
			cg_array_write("Reynolds", CG_RealDouble,1,&one, &refQ.Re);
			cg_goto(f,iBase,"ReferenceState_t",1, "DataArray_t",++count,NULL);
			cg_dataclass_write( CG_NondimensionalParameter );
		}

		// Temperature
		if( ! isnan(refQ.Tinf) )
		{
			cg_goto(f,iBase,"ReferenceState_t",1, NULL);
			cg_array_write("Temperature", CG_RealDouble,1,&one, &refQ.Tinf);
			cg_goto(f,iBase,"ReferenceState_t",1, "DataArray_t",++count,NULL);
			cg_dataclass_write( CG_Dimensional );
		}
	}

	//
	// Flow equation set
	//
	cg_goto(f,iBase, NULL);
	if( cg_equationset_write(G_Domain.nDim) == CG_OK )
	{
		const cgsize_t one = 1;
		std::map<std::string, double>::const_iterator iterPrm;

		// Gas model
		cg_goto(f,iBase,"FlowEquationSet_t",1, NULL);
		cg_model_write("GasModel_t", CG_CaloricallyPerfect );

		// TODO: Chemical kinetics
		iterPrm = G_Domain.mapCasePrms_real.find("gamma");
		if( iterPrm != G_Domain.mapCasePrms_real.end() )
		{
			// Cp/Cv
			const double gamma = iterPrm->second;
			cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, NULL);
			cg_array_write("SpecificHeatRatio", CG_RealDouble,1,&one, &gamma);
			cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, "DataArray_t",1, NULL);
			cg_dataclass_write( CG_NondimensionalParameter );

			// Ideal Gas Constant R = 1/(gamma*M^2) in nondimensional case
			{
				const double& M = G_Domain.phys->ref_state.M;
				const double R = 1./(gamma * M*M);

				cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, NULL);
				cg_array_write("IdealGasConstant", CG_RealDouble,1,&one, &R);
				cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, "DataArray_t",2, NULL);
				cg_dataclass_write( CG_NondimensionalParameter );
			}
		}


		// Thermal Conductivity Model
		cg_goto(f,iBase,"FlowEquationSet_t",1, NULL);
		cg_model_write("ThermalConductivityModel_t", CG_ConstantPrandtl);

		iterPrm = G_Domain.mapCasePrms_real.find("Pr");
		if( iterPrm != G_Domain.mapCasePrms_real.end() )
		{
			const double Pr = iterPrm->second;
			cg_goto(f,iBase,"FlowEquationSet_t",1, "ThermalConductivityModel_t",1, NULL);
			cg_array_write("Prandtl", CG_RealDouble,1,&one, &Pr);
			cg_goto(f,iBase,"FlowEquationSet_t",1, "ThermalConductivityModel_t",1, "DataArray_t",1, NULL);
			cg_dataclass_write( CG_NondimensionalParameter );
		}
	}


	// Solution time
	{
		cg_biter_write(f,iBase, "TimeIterValues", 1/*time steps count*/);
		cg_goto(f,iBase,"BaseIterativeData_t",1, NULL);

		const hsflow::TSCDtab& dscr_tab = G_Plugins.get_discret_caps().scheme->tab;
		const double& dt = dscr_tab.htau;
		const double t = dscr_tab.tau - time_layer * dt;

		const cgsize_t len = 1;
		cg_array_write("TimeValues", CG_RealDouble,1,&len, &t);
		cg_array_write("TimeStepValues"/*non-standart*/, CG_RealDouble,1,&len, &dt);
	}

	return true;
}
//-----------------------------------------------------------------------------


/**
 * Generate unique face name for CGNS
 *
 * @param[in] zone_idx - zone index, 0-based
 * @param[in] face_pos - face position
 *
 * @return face name appropriate for CGNS
 */
static std::string cgns_face_name(const int zone_idx, const TZoneFacePos face_pos)
{
	const TZone& zne = G_Domain.Zones[zone_idx];
	const int face_idx = static_cast<int>(face_pos);

	std::string name = hs_string_format("%s-%d", zne.szName, face_idx+1);
	if( name.length() > 32 )
		name = hs_string_format("%s-z%d-bc%d", zne.Faces[face_idx].szBC,
			G_Domain.map_iZne2cgID[zone_idx], face_idx+1 );

	return name;
}
//-----------------------------------------------------------------------------
