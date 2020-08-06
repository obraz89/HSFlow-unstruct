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
bool initField_struct()
{

	//
	// Memory allocation for field data
	//

	// Space to store field data of zones local to current MPI rank
	for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
	{
		TZone& zne = G_Domain.Zones[b];
		const int N = G_Domain.nu * zne.nodes_count();

		hsLogMessage("struct Zone#%d NNodes: %d", b, N);

		// Current (n+1) time layer
		zne.U = new double[N];

	}

	//hsLogDebug( "Field allocated. Memory: %.2f MiB", getCurrentRSS_MiB() );

	double time = -1;

	//
	// Fill-in uniform field data
	//
	{
		const double Uinf[] = {0,0,0,0,0};

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

	return true;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


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
bool saveField_struct(const std::string& fileName, const std::string& gridFileName,
			   const short time_layer, bool isDouble)
{

	int f = -1, fGrid = -1;  // file descriptors for field & grid
	std::string fn_grid;  // grid file name relative to the field
	int iBase = -1,  iBaseGrid = -1;  // CGNS base ids

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
	//if( ! writeMetaInfoToCGNS(f, iBase, time_layer) )
	//	break;

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

	MPI_Bcast(&ok, 1, MPI_SHORT, 0/*root*/, MPI_COMM_WORLD);
	if( ! ok )
		return false;

	// Inform all ranks if we need grid coordinates or not
	MPI_Bcast(&fGrid, 1, MPI_INT, 0/*root*/, MPI_COMM_WORLD);

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
						MPI_Ssend(Vals.data(), nxyz0, Vals.type().mpi, 0/*root*/, mpiTag, MPI_COMM_WORLD);
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
							rankSrc, mpiTag, MPI_COMM_WORLD,
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
		// removed


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
					MPI_Ssend(Vals.data(), nxyz0, Vals.type().mpi, 0/*root*/, mpiTag, MPI_COMM_WORLD);
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
						rankSrc, mpiTag, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE  // don't use NULL as status, MPI_Ssend may get stuck (i.e. in Intel MPI 5.0.1)
					);
				}

				int iField = -1;
				const char* name = g_VecFuncNames[fun].c_str();
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

//-----------------------------------------------------------------------------


