///////////////////////////////////////////////////////////////////////////////
// Name:        grid-load-cgns.cpp
// Purpose:     Load grid from file in CGNS format, 2D & 3D case
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#include <array>
#include "cgnslib_wrapper.h"

#include "common_data_struct.h"
#include "common_procs_struct.h"
#include "settings.h"

//#include "grid.h"
#include "grid-load-cgns-struct.h"

#include "logging.h"

#include "mpi.h"

#define PETSC_COMM_WORLD MPI_COMM_WORLD

#if defined(_MSC_VER)
	#define strcasecmp  _stricmp
#endif
//-----------------------------------------------------------------------------

//
// Forward declarations
//
static bool parseConnectivity( TcgnsContext& ctx );  // 1-to-1 connectivity
static bool parseBCs( TcgnsContext& ctx );     // boundary conditions
static bool parseVCs( TcgnsContext& ctx );     // volume conditions (frozen zones)

static bool loadGridCoords(TcgnsContext& ctx);
//-----------------------------------------------------------------------------

// MPI_Datatype for cgsize_t from CGNS
#if ( CG_SIZEOF_SIZE == 32 )
	#define MPI_CGSIZE  MPI_INT32_T
#else
	#define MPI_CGSIZE  MPI_INT64_T
#endif
//-----------------------------------------------------------------------------


bool doLoadGrid_cgns( const std::string& gridFN )
{
	const bool is3D = ( G_Domain.nDim == 3 );

	TcgnsContext ctx;
	if( G_State.mpiRank == 0 ) do
	{
		G_Domain.nZones = 0;

		if( cg_open( gridFN.c_str(),CG_MODE_READ, &ctx.fileID ) != CG_OK )
		{
			hsLogError( "Can't open grid file '%s' for reading (%s)",
				gridFN.c_str(), cg_get_error()   );
			break;
		}

		ctx.iBase = 1; // assume only one base

		//
		// Space dimensions
		//
		char szName[CG_MAX_NAME_LENGTH];
		int dimCell = 0, dimPhys = 0;
		cg_base_read(ctx.fileID,ctx.iBase,  szName,&dimCell,&dimPhys);
		if( dimCell != G_Domain.nDim )
		{
			hsLogError( "The grid is not for %dD problems", G_Domain.nDim );
			break;
		}

		//
		// Number of zones
		//
		cg_nzones(ctx.fileID,ctx.iBase, &G_Domain.nZones);
		if( G_Domain.nZones < 1 )
		{
			hsLogError( "Domain zones are not found" );
			break;
		}
	} while(false);

	MPI_Bcast( &(G_Domain.nZones), 1, MPI_INT, 0/*root*/, PETSC_COMM_WORLD );
	if( G_Domain.nZones < 1 )
		return false;

	// Allocate memory for whole computational domain
	G_Domain.Zones = new TZone[ G_Domain.nZones ];

	// Temporary zones data used on CGNS file parsing
	ctx.cgZones = new TcgnsZone[ G_Domain.nZones ];

	// Distribute grid zones through MPI ranks (processes)
	if( ! assignZonesToProcs() )
		return false;


	//
	// Get zones sizes and names on root rank
	//
	for( int iZne = 0;  iZne < G_Domain.nZones;  ++iZne )
	{
		const int& cgZneID = G_Domain.map_iZne2cgID[iZne];

		char isOk = 0;
		if( G_State.mpiRank == 0 )
		{
			CG_ZoneType_t type;  cg_zone_type(ctx.fileID,ctx.iBase,cgZneID, &type);
			isOk = (type == CG_Structured) ? 1 : 0;
		}
		MPI_Bcast(&isOk, 1,MPI_CHAR, 0, PETSC_COMM_WORLD);
		if( ! isOk )
		{
			hsLogError( "Only structured grids are supported" );
			return false;
		}

		TZone& zne = G_Domain.Zones[iZne];

		// Get zone name & size
		int nxyz[3] = { 0, 0, 0 };
		if( G_State.mpiRank == 0 )
		{
			cgsize_t isize[9]; // 3D: NVertex{I,J,K}, NCell{I,J,K}, NBoundVertex{I,J,K}
							   // 2D: NVertex{I,J},   NCell{I,J}, NBoundVertex{I,J}
			if( cg_zone_read(ctx.fileID,ctx.iBase, cgZneID,  zne.szName,isize) == CG_OK )
			{
				nxyz[0] = isize[0];
				nxyz[1] = isize[1];
				nxyz[2] = is3D ? isize[2] : 1;
			}
		}
		MPI_Bcast(zne.szName, 40,MPI_CHAR, 0, PETSC_COMM_WORLD);
		if( ! zne.szName[0] ) {
			hsLogError( "Can't read zone #%d info", cgZneID );
			return false;
		}

		MPI_Bcast(nxyz, 3,MPI_INT, 0, PETSC_COMM_WORLD);
		zne.nx = nxyz[0];
		zne.ny = nxyz[1];
		zne.nz = nxyz[2];

		// Indices of real (not ghost) nodes
		zne.is = 1;  zne.ie = zne.nx;
		zne.js = 1;  zne.je = zne.ny;
		zne.ks = 1;  zne.ke = zne.nz;

		ctx.map_ZneName2Idx[zne.szName] = iZne;
	}


	// Volume conditions info (frozen zones)
	parseVCs(ctx);

	// Connectivity info
	// Updates {zne,cgZne}.{is,ie,js,je,ks,ke}, zne.{nx,ny,nz}
	if( ! parseConnectivity(ctx) )
		return false;

	// Boundary conditions info
	if( ! parseBCs(ctx) )
		return false;

	// Read grid coordinates
	if( ! loadGridCoords(ctx) )
		return false;

	cg_close( ctx.fileID );

	return true;
}
//-----------------------------------------------------------------------------


/**
 * Read from CGNS file & parse a 1to1 connectivity patch
 *
 * @param[in] ctx - context of the opened CGNS file
 * @param[in] iZne - 0-based internal index of the zone to read patch from
 * @param[in] cgPatchID - patch ID to process
 *
 * @param[out] cgFacePatch - parsed data
 *
 * @return true if succeded, false otherwise
 */
static bool parse_1to1_connectivity_patch( const TcgnsContext& ctx,
	const int iZne, const int cgPatchID, TcgnsZone::TFacePatch& cgFacePatch )
{
	const int& cgZneID = G_Domain.map_iZne2cgID[iZne];

	// Patch names at current & donor zones
	TpakArrays<char,2, CG_MAX_NAME_LENGTH>  pkNames;
	char *szName = pkNames[0], *szDonor = pkNames[1];

	TpakArrays<cgsize_t,2, 6> pkIndices;
	cgsize_t *idxRng = pkIndices[0],  *idxRngDnr = pkIndices[1];

	int iTsh[3]; // short-hand notation of transform matrix

	if( G_State.mpiRank == 0 )
	{
		szName[0] = 0; // error indicator
		if( cg_1to1_read( ctx.fileID,ctx.iBase,cgZneID,cgPatchID,
			              szName, szDonor, idxRng, idxRngDnr, iTsh ) != CG_OK )
		{
			hsLogError( "Can't read connection patch '%s'(#%d) of zone '%s'(#%d)",
				szName, cgPatchID, G_Domain.Zones[iZne].szName, cgZneID   );
		}
	}

	MPI_Bcast(pkNames.data(), pkNames.size(), MPI_CHAR, 0, PETSC_COMM_WORLD);
	if( ! szName[0] )
		return false;

	MPI_Bcast(pkIndices.data(), pkIndices.size(), MPI_CGSIZE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(iTsh, 3, MPI_INT, 0, PETSC_COMM_WORLD);

	//---
	const bool is3D = (G_Domain.nDim == 3);

	int m = 0;
	const cgsize_t
		&is = idxRng[m++],  &js = idxRng[m++],  &ks = is3D ? idxRng[m++] : 1,
		&ie = idxRng[m++],  &je = idxRng[m++],  &ke = is3D ? idxRng[m++] : 1;

	m = 0;
	const cgsize_t
		&ids = idxRngDnr[m++],  &jds = idxRngDnr[m++],   &kds = is3D ? idxRngDnr[m++] : 1,
		&ide = idxRngDnr[m++],  &jde = idxRngDnr[m++],   &kde = is3D ? idxRngDnr[m++] : 1;

	cgFacePatch.is0 = is;   cgFacePatch.ie0 = ie;
	cgFacePatch.js0 = js;   cgFacePatch.je0 = je;
	cgFacePatch.ks0 = ks;   cgFacePatch.ke0 = ke;

	//
	// Detect face, the current connection patch belongs to
	//
		 if( is==ie )  cgFacePatch.posFace = (is==1) ? faceXmin : faceXmax;
	else if( js==je )  cgFacePatch.posFace = (js==1) ? faceYmin : faceYmax;
	else if( ks==ke )  cgFacePatch.posFace = (ks==1) ? faceZmin : faceZmax;

	//
	// Detect donor zone
	//
	const auto& iterZneDnr = ctx.map_ZneName2Idx.find(szDonor);
	if( iterZneDnr == ctx.map_ZneName2Idx.end() )
	{
		hsLogError(
			"Connection patch '%s'#%d of zone '%s'#%d abuts to nonexistent zone '%s'",
			szName, cgPatchID, G_Domain.Zones[iZne].szName, cgZneID,  szDonor );
		return false;
	}
	cgFacePatch.iZneDnr = iterZneDnr->second;

	//
	// Detect face, the donor connection patch belongs to.
	// Check if connectivity info is correct
	//
	const TZone& zneDnr = G_Domain.Zones[cgFacePatch.iZneDnr];
	if( ids == ide )
	{
		if( ids == 1 )
			cgFacePatch.dnrPosFace = faceXmin;
		else if( ids == (zneDnr.ie - zneDnr.is + 1) ) // original count
			cgFacePatch.dnrPosFace = faceXmax;
	}
	else if( jds == jde )
	{
		if( jds == 1 )
			cgFacePatch.dnrPosFace = faceYmin;
		if( jds == (zneDnr.je - zneDnr.js + 1) )
			cgFacePatch.dnrPosFace = faceYmax;
	}
	else if( kds == kde )
	{
		assert(is3D);
		if( kds == 1 )
			cgFacePatch.dnrPosFace = faceZmin;
		else if( kds == (zneDnr.ke - zneDnr.ks + 1) )
			cgFacePatch.dnrPosFace = faceZmax;
	}

	if( cgFacePatch.dnrPosFace == faceNone ) {
		hsLogError(
			"Connection patch '%s'#%d of zone '%s'#%d doesn't abut to any face",
			szName, cgPatchID, G_Domain.Zones[iZne].szName, cgZneID   );
		return false;
	}
	cgFacePatch.dnr_is0 = ids;
	cgFacePatch.dnr_js0 = jds;
	cgFacePatch.dnr_ks0 = kds;


	// Transform matrix between node indexes
	// of the current and donor zones
	//
	// CGNS SIDS documentation,
	// section 8.2 1-to-1 Interface Connectivity Structure Definition
	//
	// iTsh = [a, b, c] =>
	//	   | sgn(a)del(a,1) sgn(b)del(b,1) sgn(c)del(c,1) |
	// T = | sgn(a)del(a,2) sgn(b)del(b,2) sgn(c)del(c,2) |,
	//	   | sgn(a)del(a,3) sgn(b)del(b,3) sgn(c)del(c,3) |
	auto sgn = [](int x){ return (x<0) ? -1 : +1; };
	auto del = [](int x, int y){ return (abs(x)==abs(y)) ? +1 : 0; };

	int* mT = cgFacePatch.matTrans;
	mT[0] =        sgn(iTsh[0]) * del(iTsh[0],1);
	mT[1] =        sgn(iTsh[1]) * del(iTsh[1],1);
	mT[2] = is3D ? sgn(iTsh[2]) * del(iTsh[2],1) : 0;

	mT[3] =        sgn(iTsh[0]) * del(iTsh[0],2);
	mT[4] =        sgn(iTsh[1]) * del(iTsh[1],2);
	mT[5] = is3D ? sgn(iTsh[2]) * del(iTsh[2],2) : 0;

	mT[6] =        sgn(iTsh[0]) * del(iTsh[0],3);
	mT[7] =        sgn(iTsh[1]) * del(iTsh[1],3);
	mT[8] = is3D ? sgn(iTsh[2]) * del(iTsh[2],3) : 0;


	//
	// Periodic BC
	//
	TZoneFace& face = G_Domain.Zones[iZne].Faces[cgFacePatch.posFace];
	if( face.transform )
		return true;  // already filled by other patch

	TpakArrays<float,3, 3> pkTransform;
	float *C = pkTransform[0], // rotation center
		  *A = pkTransform[1], // rotation angle
		  *T = pkTransform[2]; // translation vector

	if( G_State.mpiRank == 0 ) {
		C[0] = NAN;  // error indicator
		cg_1to1_periodic_read(ctx.fileID, ctx.iBase, cgZneID, cgPatchID, C, A, T);
	}
	MPI_Bcast(pkTransform.data(), pkTransform.size(), MPI_FLOAT, 0, PETSC_COMM_WORLD );
	if( isnan(C[0]) )
		return true;

	// Reversed transform (from donor to the current face)
	A[0] = -A[0];   A[1] = -A[1];   A[2] = -A[2];
	T[0] = -T[0];   T[1] = -T[1];   T[2] = -T[2];

	face.transform = new TZoneFace::Ttransform;
	copy_n(C, 3, face.transform->TransformCenter);
	copy_n(T, 3, face.transform->Translation);

	// Rotation matrix
	double* R = face.transform->TransformMatrix;
	if( ! is3D )
	{
		const double g = A[0];  // gamma
		R[0] = cos(g);   R[1] = -sin(g);
		R[2] = sin(g);   R[3] = cos(g);

		R[4] = R[5] = R[6] = R[7] = R[8] = NAN;
	}
	else
	{
		const double a = A[0],  b = A[1],  g = A[2];  // alpha, beta, gamma
		R[0] = cos(b)*cos(g);   R[1] = cos(g)*sin(a)*sin(b) - cos(a)*sin(g);   R[2] = cos(a)*cos(g)*sin(b) + sin(a)*sin(g);
		R[3] = cos(b)*sin(g);   R[4] = cos(a)*cos(g) - sin(a)*sin(b)*sin(g);   R[5] = -cos(g)*sin(a) + cos(a)*sin(b)*sin(g);
		R[6] = -sin(b);         R[7] = cos(b)*sin(a);                          R[8] = cos(a)*cos(b);
	}

	//
	// face.transform->{is,ie,js,je,ks,ke} are to be assigned
	//

	return true;
}
//-----------------------------------------------------------------------------


/**
 *  Read connectivity data from the opened CGNS file on root MPI rank
 *  and initialize ghost indices
 *
 *  @param[in] ctx - context of the opened CGNS file
 *
 *  @return true if succeeded, false if failed
**/
static bool parseConnectivity( TcgnsContext& ctx )
{
	const bool is3D = (G_Domain.nDim == 3);

	// Number of ghost nodes (half of stencil size)
	const int ghostI = 2, ghostJ = 2, ghostK = 2;
	// TODO:
	//G_Plugins.get_discret_caps().getParI(11/*nsx*/, ghostI);  ghostI /= 2;
	//G_Plugins.get_discret_caps().getParI(24/*nsy*/, ghostJ);  ghostJ /= 2;

//
// (1) Read connectivity data from CGNS file
//
for( int iZne = 0;  iZne < G_Domain.nZones;  ++iZne )
{
	const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
	TZone& zne = G_Domain.Zones[iZne];
	TcgnsZone& cgZne = ctx.cgZones[iZne];

	if( G_State.mpiRank == 0 )
		cg_n1to1(ctx.fileID,ctx.iBase,cgZneID, &cgZne.nPatches);
	MPI_Bcast( &cgZne.nPatches, 1, MPI_INT, 0, PETSC_COMM_WORLD );

	if( cgZne.nPatches < 1 && G_Domain.nZones > 1 ) {
		hsLogError("Missing inter-zone 1-to-1 connectivity in '%s'#%d", zne.szName, cgZneID);
		return false;
	}
	cgZne.Patches = new TcgnsZone::TFacePatch[cgZne.nPatches];

	// Loop through connectivity patches
	for( int cgPatchId = 1; cgPatchId <= cgZne.nPatches; ++cgPatchId )
	{
		TcgnsZone::TFacePatch& cgFacePatch = cgZne.Patches[cgPatchId - 1];
		if( ! parse_1to1_connectivity_patch(ctx,iZne, cgPatchId, cgFacePatch) )
			return false;

		//
		// Increase the zone size by ghosts count
		//
		auto& face_type = cgZne.face_types[cgFacePatch.posFace];
		if( face_type == TcgnsZone::FaceType::unknown )
		{
			face_type = TcgnsZone::FaceType::abutted;

			if( ! zne.isFrozen )
			switch( cgFacePatch.posFace )
			{
			case faceXmin:  zne.nx += ghostI;  zne.is += ghostI;  zne.ie += ghostI;  break;
			case faceXmax:  zne.nx += ghostI;                                        break;

			case faceYmin:  zne.ny += ghostJ;  zne.js += ghostJ;  zne.je += ghostJ;  break;
			case faceYmax:  zne.ny += ghostJ;                                        break;

			case faceZmin:  zne.nz += ghostK;  zne.ks += ghostK;  zne.ke += ghostK;  break;
			case faceZmax:  zne.nz += ghostK;                                        break;
			}
		}

	} // loop connection patches
} // loop through zones


//
// (2) Assign owners for abutted faces
//

//
// (2.1) Mark faces that should be skipped and processed by abutted zone
//
for( int b = 0; b < G_Domain.nZones; ++b )
{
	TZone& zne = G_Domain.Zones[b];
	const TcgnsZone& cgZne = ctx.cgZones[b];

	// Don't skip any grid layer from a frozen zone,
	// it has no ghost nodes to substitute.
	if( zne.isFrozen )
		continue;

	bool face_processed[6] = {};  // inits all by false
	for( int p = 0; p < cgZne.nPatches; ++p ) // face patches
	{
		const TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[p];
		const TZone& zneDnr = G_Domain.Zones[ cgPatch.iZneDnr ];
		const TcgnsZone& cgZneDnr = ctx.cgZones[cgPatch.iZneDnr];

		if( face_processed[cgPatch.posFace] )
			continue;

		// Step from frozen zone
		if( zneDnr.isFrozen )
		{
			zne.Faces[cgPatch.posFace].isSkipped = true;
			continue;
		}

		//
		// Skip current face if donor face adjacent to more boundaries
		//
		auto adjFaces = [](const TZoneFacePos& f) -> std::array<TZoneFacePos,4>
		{
			switch( f ) {
			case faceXmin:  case faceXmax:
				return {{faceYmin, faceYmax, faceZmin, faceZmax}};
			case faceYmin:   case faceYmax:
				return {{faceXmin, faceXmax, faceZmin, faceZmax}};
			case faceZmin:   case faceZmax:
				return {{faceXmin, faceXmax, faceYmin, faceYmax}};
			default:
				return {{faceNone, faceNone, faceNone, faceNone}};
			}
		};
		auto af  = adjFaces(cgPatch.posFace);    // adjacent faces
		auto afD = adjFaces(cgPatch.dnrPosFace); // donor adjacent faces

		// Count of adjacent faces with BC
		int BCcount = 0, BCcountD = 0;
		for( int f = 0;  f < 2* (G_Domain.nDim - 1);  ++f )
		{
			if( cgZne.face_types[af[f]] != TcgnsZone::FaceType::abutted )
				++BCcount;

			if( cgZneDnr.face_types[afD[f]] != TcgnsZone::FaceType::abutted )
				++BCcountD;
		}

		if( BCcountD > BCcount )
			zne.Faces[cgPatch.posFace].isSkipped = true;

		face_processed[cgPatch.posFace] = true;
	} // loop face patches
} // loop zones


#if 0
//
// (2.2) Skip faces owned by both abutted zones
//
// FIXME: Can we skip some corner nodes entirely this way ???
//
for( int zi = 0; zi < G_Domain.nZones; ++zi )
{
	TcgnsZone& cgZne = ctx.cgZones[zi];
	for( int p = 0; p < cgZne.nPatches; ++p ) // abutted face patches
	{
		TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[p];

		TZoneFace& face = G_Domain.Zones[zi].Faces[cgPatch.posFace];
		TZoneFace& faceDnr = G_Domain.Zones[cgPatch.iZneDnr].Faces[cgPatch.dnrPosFace];

		if( ! face.isSkipped && ! faceDnr.isSkipped )
			faceDnr.isSkipped = true;
	}
}
#endif

//
// (3) Update domain zones dimensions and indices using abutted faces info
//
for( int b = 0; b < G_Domain.nZones; ++b )
{
	TZone& zne = G_Domain.Zones[b];
	TcgnsZone& cgZne = ctx.cgZones[b];

	cgZne.is1 = zne.is;   cgZne.ie1 = zne.ie;
	cgZne.js1 = zne.js;   cgZne.je1 = zne.je;
	cgZne.ks1 = zne.ks;   cgZne.ke1 = zne.ke;

	for( int f = 0; f < 2*G_Domain.nDim; ++f ) // faces
	{
		if( cgZne.face_types[f] != TcgnsZone::FaceType::abutted )
			continue;

		TZoneFace& face = zne.Faces[f];
		if( face.isSkipped )
		{
			// Face is owned by donor zone, skip it here -> cut nodes layer
			switch( f )
			{
			case faceXmin:
				zne.nx -= 1;
				cgZne.ie1 = (zne.ie -= 1);
				cgZne.is1 -= 1;
				break;
			case faceXmax:
				zne.nx -= 1;
				zne.ie -= 1;
				break;
			case faceYmin:
				zne.ny -= 1;
				cgZne.je1 = (zne.je -= 1);
				cgZne.js1 -= 1;
				break;
			case faceYmax:
				zne.ny -= 1;
				zne.je -= 1;
				break;
			case faceZmin:
				zne.nz -= 1;
				cgZne.ke1 = (zne.ke -= 1);
				cgZne.ks1 -= 1;
				break;
			case faceZmax:
				zne.nz -= 1;
				zne.ke -= 1;
				break;
			}
		}

		// Indices range of nodes to transform due to periodicity
		if( face.transform )
		{
			TZoneFace::Ttransform& t = *(face.transform);
			t.is = 1;   t.ie = zne.nx;
			t.js = 1;   t.je = zne.ny;
			t.ks = 1;   t.ke = zne.nz;

			switch( f )
			{
			case faceXmin:  t.is = 1;          t.ie = zne.is-1;   break;
			case faceXmax:  t.is = zne.ie+1;   t.ie = zne.nx;     break;
			case faceYmin:  t.js = 1;          t.je = zne.js-1;   break;
			case faceYmax:  t.js = zne.je+1;   t.je = zne.ny;     break;
			case faceZmin:  t.ks = 1;          t.ke = zne.ks-1;   break;
			case faceZmax:  t.ks = zne.ke+1;   t.ke = zne.nz;     break;
			}
		}

	} // loop faces
} // loop zones



//
// (4) Global indices of real nodes
//
int gIdx = 0;
for( int b = 0; b < G_Domain.nZones; ++b )
{
	TZone& zne = G_Domain.Zones[b];
	TcgnsZone& cgZne = ctx.cgZones[b];
	cgZne.globIndices = new TZoneGlobIndices(zne);

	zne.nGlobStart = gIdx;

	// Real nodes count (without ghosts)
	const int
		nx0 = zne.ie - zne.is + 1,
		ny0 = zne.je - zne.js + 1,
		nz0 = zne.ke - zne.ks + 1;

	gIdx += nx0*ny0*nz0;
}


//
// (5) Parse connectivity data and fill in ghost nodes indices
//
for( int stage = 0; stage < G_Domain.nDim; ++stage ){
for( int b = 0; b < G_Domain.nZones; ++b )
{
	TZone& zne = G_Domain.Zones[b];
	TcgnsZone& cgZne = ctx.cgZones[b];
	TZoneGlobIndices& zneGlobInd = *(cgZne.globIndices);

	if( zne.isFrozen )  continue;

	// Loop through face patches
	for( int p = 0; p < cgZne.nPatches; ++p )
	{
		const TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[p];

		// Start indices of the original patch (assuming no layers were skipped) in working numbering (with ghosts)
		const int ips = cgPatch.is0 + cgZne.is1 - 1;
		const int jps = cgPatch.js0 + cgZne.js1 - 1;
		const int kps = cgPatch.ks0 + cgZne.ks1 - 1;

		//
		// Start & End indices of ghost layer
		int igs, ige,   jgs, jge,   kgs, kge;

		// Indices range of the patch in working numbering assuming no layers were skipped
		igs = cgPatch.is0 + cgZne.is1 - 1;
		ige = cgPatch.ie0 + cgZne.is1 - 1;
		jgs = cgPatch.js0 + cgZne.js1 - 1;
		jge = cgPatch.je0 + cgZne.js1 - 1;
		kgs = cgPatch.ks0 + cgZne.ks1 - 1;
		kge = cgPatch.ke0 + cgZne.ks1 - 1;

		if( stage > 0 )
		{
			// Extend patch nodes by ghosts if patch edge coinside with the face edge
			// NB: patch index range may be unsorted
			int nn = cgZne.ie1 - cgZne.is1 + 1;
			     if( cgPatch.is0 == 1  )  igs = 1;
			else if( cgPatch.is0 == nn )  igs = zne.nx;
			     if( cgPatch.ie0 == 1  )  ige = 1;
			else if( cgPatch.ie0 == nn )  ige = zne.nx;

			nn = cgZne.je1 - cgZne.js1 + 1;
			     if( cgPatch.js0 == 1  )  jgs = 1;
			else if( cgPatch.js0 == nn )  jgs = zne.ny;
			     if( cgPatch.je0 == 1  )  jge = 1;
			else if( cgPatch.je0 == nn )  jge = zne.ny;

			nn = cgZne.ke1 - cgZne.ks1 + 1;
			     if( cgPatch.ks0 == 1  )  kgs = 1;
			else if( cgPatch.ks0 == nn )  kgs = zne.nz;
			     if( cgPatch.ke0 == 1  )  kge = 1;
			else if( cgPatch.ke0 == nn )  kge = zne.nz;
		}

		// Nodes of the ghost layer extruded from the patch
		switch( cgPatch.posFace )
		{
		case faceXmin:  igs = 1;          ige = zne.is-1;   break;
		case faceXmax:  igs = zne.ie+1;   ige = zne.nx;     break;
		case faceYmin:  jgs = 1;          jge = zne.js-1;   break;
		case faceYmax:  jgs = zne.je+1;   jge = zne.ny;     break;
		case faceZmin:  kgs = 1;          kge = zne.ks-1;   break;
		case faceZmax:  kgs = zne.ke+1;   kge = zne.nz;     break;
		}

		const TZone& zneDnr = G_Domain.Zones[ cgPatch.iZneDnr ];
		const TcgnsZone& cgZneDnr = ctx.cgZones[ cgPatch.iZneDnr ];
		TZoneGlobIndices& zneDnrGlobInd = *(cgZneDnr.globIndices);

		// Start indices of the original donor patch (assuming no layers were skipped) in working numbering
		const int ids = cgPatch.dnr_is0 + cgZneDnr.is1 - 1;
		const int jds = cgPatch.dnr_js0 + cgZneDnr.js1 - 1;
		const int kds = cgPatch.dnr_ks0 + cgZneDnr.ks1 - 1;

		// Swap range limits if not ascending
		auto swap = [](int& s, int& e){  if(s > e){ int s0 = s;  s = e;  e = s0; } };
		swap(igs,  ige);   swap(jgs,  jge);   if(is3D) swap(kgs,  kge);

		for( int k = kgs; k <= kge; ++k ){
		for( int j = jgs; j <= jge; ++j ){
		for( int i = igs; i <= ige; ++i )
		{
			// Donor zone indices
			int id = (i-ips)*cgPatch.matTrans[0] + (j-jps)*cgPatch.matTrans[1] + (k-kps)*cgPatch.matTrans[2] + ids;
			int jd = (i-ips)*cgPatch.matTrans[3] + (j-jps)*cgPatch.matTrans[4] + (k-kps)*cgPatch.matTrans[5] + jds;
			int kd = (i-ips)*cgPatch.matTrans[6] + (j-jps)*cgPatch.matTrans[7] + (k-kps)*cgPatch.matTrans[8] + kds;

			if( 0 == stage )
			{
				if( id < cgZneDnr.is1 || id > cgZneDnr.ie1 ||
					jd < cgZneDnr.js1 || jd > cgZneDnr.je1 ||
					kd < cgZneDnr.ks1 || kd > cgZneDnr.ke1
					)
				{
					hsLogError(
						"Ghost node (%d,%d,%d) of zone '%s' tried to be mapped "
						"to nonexistent node (%d,%d,%d) of zone '%s'. Check indices orientation!",
						i,  j,  k,  zne.szName,
						id, jd, kd, zneDnr.szName
					);
					return false;
				}

				if( id < zneDnr.is || id > zneDnr.ie ||
					jd < zneDnr.js || jd > zneDnr.je ||
				    kd < zneDnr.ks || kd > zneDnr.ke    )
				{
					continue;    // process on next stage
				}

				zneGlobInd(i,j,k) = zneDnr.globRealInd(id,jd,kd);
			}
			else // resolve cross references (nodes around corners)
			{
				if( zneGlobInd(i,j,k) != -1 )  continue;

				if( id < 1 || id > zneDnr.nx ||
					jd < 1 || jd > zneDnr.ny ||
				    kd < 1 || kd > zneDnr.nz    )
				{
					// ghost data unavailable
					continue;
				}

				zneGlobInd(i,j,k) = zneDnrGlobInd(id,jd,kd);
			}
		}}} // for i,j,k
	} // loop face patches
}} // loop zones, stages


//
// Fill-in globIndices[] of zones in current MPI-rank
//
for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
{
	TZone& zne = G_Domain.Zones[b];
	TcgnsZone& cgZne = ctx.cgZones[b];
	TZoneGlobIndices& zneGlobInd = *(cgZne.globIndices);

	try {
		zne.globIndices = new int[ zne.nodes_count() ];
	}
	catch( std::bad_alloc& e )
	{
		printf(
			"[%d] %s(): Memory allocation for TZone::globIndices[] in zone #%d/%d FAILED. Asked for %.3f MB.\n",
			G_State.mpiRank, __FUNCTION__,
			b+1, G_Domain.nZones,
			zne.nodes_count() * (float)sizeof(int) / (1024.*1024.)
		);
		return false;
	}

	// Fill-in global indices of real nodes
	for( int k=1; k<=zne.nz; ++k ){
	for( int j=1; j<=zne.ny; ++j ){
	for( int i=1; i<=zne.nx; ++i )
	{
		zne.globInd(i,j,k) = zneGlobInd(i,j,k);
	}}}

#if ! defined(NDEBUG) && 0
	debug_connectivity(b);
#endif
}

	return true;
}
//-----------------------------------------------------------------------------


/**
 * Parse CGNS boundary condition info
 *
 * @pre All abutted faces should be processed already & ctx.cgZones[] filled
 *
 * @param ctx - info concerning the CGNS file
 * @return true if succeeded
 */
static bool parseBCs( TcgnsContext& ctx )
{
	const bool is3D = (G_Domain.nDim == 3);

for( int iZne = 0;  iZne < G_Domain.nZones;  ++iZne )
{
	const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
	TZone& zne = G_Domain.Zones[iZne];
	TcgnsZone& cgZne = ctx.cgZones[iZne];

	// Original zone size (not accounting added ghosts & skipped layers)
	const int
		nx0 = cgZne.ie1 - cgZne.is1 + 1,
		ny0 = cgZne.je1 - cgZne.js1 + 1,
		nz0 = cgZne.ke1 - cgZne.ks1 + 1;

	// Number of BCs in the Zone
	int nBCs = 0;
	if( G_State.mpiRank == 0 )   cg_nbocos(ctx.fileID,ctx.iBase,cgZneID, &nBCs);
	MPI_Bcast(&nBCs, 1,MPI_INT, 0, PETSC_COMM_WORLD);

	for( int iBC = 1; iBC <= nBCs; ++iBC )
	{
		char szPatchName[CG_MAX_NAME_LENGTH];
		if( G_State.mpiRank == 0 )
		{
			CG_BCType_t iBCtype;

			CG_PointSetType_t pntSetType;
			cgsize_t nPnts = 0; // number of points defining the BC region

			// Normals to the BC patch - UNUSED
			int iNorm[3]; // normal as index vector (computational coords)
			cgsize_t normListSize = 0;  CG_DataType_t normDataType; // normals in phys coords

			int nDatasets = 0; // number of datasets with additional info for the BC

			cg_boco_info( ctx.fileID, ctx.iBase, cgZneID, iBC,
				szPatchName, &iBCtype,
				&pntSetType, &nPnts,
				iNorm, &normListSize, &normDataType,
				&nDatasets
			);
			if( pntSetType != CG_PointRange && nPnts != 2 )
			{
				hsLogError(
					"Boundary condition patch '%s'(#%d) of zone '%s'(#%d) isn't defined as point range",
					szPatchName, iBC,   zne.szName, cgZneID   );
				szPatchName[0] = 0x3;  // 'end of text' code -> error indicator
			}
		}
		MPI_Bcast(szPatchName, CG_MAX_NAME_LENGTH,MPI_CHAR, 0, PETSC_COMM_WORLD);
		if( szPatchName[0] == 0x3 )
			return false;

		// Read BC patch point range
		cgsize_t idxRng[6];   // {I,J,K}min,  {I,J,K}max
		if( G_State.mpiRank == 0 )
			cg_boco_read(ctx.fileID,ctx.iBase,cgZneID, iBC, idxRng, nullptr);

		MPI_Bcast(idxRng, 6,MPI_CGSIZE, 0, PETSC_COMM_WORLD);

		int m = 0;
		const cgsize_t
			&is = idxRng[m++],  &js = idxRng[m++],  &ks = is3D ? idxRng[m++] : 1,
			&ie = idxRng[m++],  &je = idxRng[m++],  &ke = is3D ? idxRng[m++] : 1;

		// Check that BC patch cover whole face
		bool ok = true;
		ok &= (is==1 && ie==1) || (is==nx0 && ie==nx0) || (is==1 && ie==nx0);
		ok &= (js==1 && je==1) || (js==ny0 && je==ny0) || (js==1 && je==ny0);
		ok &= (ks==1 && ke==1) || (ks==nz0 && ke==nz0) || (ks==1 && ke==nz0);
		if( ! ok )
		{
			hsLogError(
				"Boundary condition patch '%s'#%d of zone '%s'#%d doesn't cover whole zone face",
				szPatchName, iBC,
				zne.szName, cgZneID
			);
			return false;
		}

		// BC family name
		char szBC[CG_MAX_NAME_LENGTH] = "";
		if( G_State.mpiRank == 0 )
		{
			cg_goto(ctx.fileID,ctx.iBase, "Zone_t",cgZneID, "ZoneBC",0, "BC_t",iBC, NULL);
			if( cg_famname_read(szBC) == CG_OK )
			{
				// Read "Fam_Descr_Name" generated by Pointwise 16.03
				if( cg_goto(ctx.fileID,ctx.iBase, szBC,0, NULL) == CG_OK )
				{
					char szDescrName[CG_MAX_NAME_LENGTH] = "";
					char* szDescrText = nullptr;
					if( cg_descriptor_read(1, szDescrName, &szDescrText) == CG_OK )
					{
						if( strcmp(szDescrName, "Fam_Descr_Name") == 0 )
						{
							strcpy(szBC, szDescrText);
							cg_free(szDescrText);
						}
					}
				}
			}
			else
			{
				strcpy(szBC, szPatchName);
			}
		}
		MPI_Bcast(szBC, CG_MAX_NAME_LENGTH,MPI_CHAR, 0, PETSC_COMM_WORLD);

		// Detect BC face
		TZoneFacePos posFace = faceNone;
		     if( is==ie )  posFace = ( is==1 ) ? faceXmin : faceXmax;
		else if( js==je )  posFace = ( js==1 ) ? faceYmin : faceYmax;
		else if( ks==ke )  posFace = ( ks==1 ) ? faceZmin : faceZmax;

		// BC name
		strcpy(zne.Faces[posFace].szBC,   szBC);

		cgZne.face_types[posFace] = TcgnsZone::FaceType::boundary;
	} // for iBC

	//
	// Fill-in missed BCs
	//
	for( int f = 0; f < G_Domain.nDim * 2; ++f )
	{
		// All known faces should be already marked as `abutted` or `boundary`
		if( cgZne.face_types[f] != TcgnsZone::FaceType::unknown )
			continue;

		TZoneFace& face = zne.Faces[f];

		const char* names[] = { "Imin", "Imax", "Jmin", "Jmax", "Kmin", "Kmax" };
		std::string bc = hs_string_format("%s-%s", zne.szName, names[f]);

		hsLogWarning(  "Zone '%s' has face w/o BC - assigning '%s'",
			zne.szName,   bc.c_str() );

		if( bc.size() > sizeof(face.szBC) ) {
			hsLogError("Too long BC name");
			return false;
		}
		strcpy(face.szBC, bc.c_str());
	}

} // for iZone

	return true;
}
//-----------------------------------------------------------------------------


/**
 * Parse CGNS volume condition (VC) data on root rank
 *
 * @param ctx info concerning the CGNS file
 * @return true if succeeded
 */
static bool parseVCs( TcgnsContext& ctx )
{
	//if( ! g_genOpts.allowFrozenZones )  return true;

	return true;
}
//-----------------------------------------------------------------------------

/**
 * Read grid coordinates on root MPI rank, and scatter to the owner ranks
 *
 * @param ctx info concerning the CGNS file
 * @return true if succeeded
 */
static bool loadGridCoords(TcgnsContext& ctx)
{
	const bool is3D = (G_Domain.nDim == 3);

	for( int iZne = 0;  iZne < G_Domain.nZones;  ++iZne )
	{
		TZone& zne = G_Domain.Zones[iZne];

		//
		// Grid step in computational (curvilinear) coordinates
		//
		// TODO: check g_COMP_SPACE_GRID_STEP
		zne.grd.dksi = 1. / (zne.nx - 1);
		zne.grd.deta = 1. / (zne.ny - 1);
		if (is3D) zne.grd.dzet = 1. / (zne.nz - 1);

		zne.grd.l_dksi = 1.0 / zne.grd.dksi;
		zne.grd.l_deta = 1.0 / zne.grd.deta;
		zne.grd.l_dzet = (is3D) ? (1.0 / zne.grd.dzet) : HUGE_VAL;

		const TcgnsZone& cgZne = ctx.cgZones[iZne];

		// Original zone size (not accounting added ghosts & skipped layers)
		const int
			nx0 = cgZne.ie1 - cgZne.is1 + 1,
			ny0 = cgZne.je1 - cgZne.js1 + 1,
			nz0 = cgZne.ke1 - cgZne.ks1 + 1;

		//
		// Load coordinates
		//
		short ok = 0;
		TpakArraysDyn<double>  pkXYZ(G_Domain.nDim, nx0*ny0*nz0);
		if( G_State.mpiRank == 0 )
		{
			const int& cgZneID = G_Domain.map_iZne2cgID[iZne];
			const char* coordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };
			cgsize_t irmin[3] = {1, 1, 1};       // indices range to read
			cgsize_t irmax[3] = {nx0, ny0, nz0}; // (all zone)

			for( int m = 0; m < G_Domain.nDim; ++m )
			{
				const char* name = coordNames[m];
				ok = ( cg_coord_read( ctx.fileID,ctx.iBase,cgZneID, name, CG_RealDouble, irmin,irmax, pkXYZ[m] ) == CG_OK );
				if( ! ok ) {
					hsLogError("Can't read %s from zone '%s'#%d (%s)",
						name,   zne.szName, cgZneID,   cg_get_error()  );
					break;
				}
			}
		}
		MPI_Bcast(&ok,1,MPI_SHORT, 0, PETSC_COMM_WORLD);
		if( ! ok )
			return false;

		//
		// Send data from root to the zone's owner
		//
		const int mpiTag = 'c'+'o'+'o'+'r'+'d' + iZne;
		if( G_State.mpiRank == 0 )
		{
			const int& rankDst = G_State.map_zone2rank[iZne];
			if( rankDst != 0 )  // don't send to myself
				MPI_Ssend( pkXYZ.data(), pkXYZ.size(), MPI_DOUBLE, rankDst, mpiTag, PETSC_COMM_WORLD );
		}
		else if( G_Domain.bs <= iZne && iZne <= G_Domain.be )
		{
			MPI_Recv( pkXYZ.data(), pkXYZ.size(), MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE );
		}

		//
		// Pack coordinates into internal structures
		//
		if( G_Domain.bs <= iZne && iZne <= G_Domain.be )
		{
			if( is3D )
				zne.grd.c3d = new Tcoord3D[zne.nodes_count()];
			else
				zne.grd.c2d = new Tcoord2D[zne.nodes_count()];

			if( is3D )
			{
				for( int k = 1; k <= zne.nz; ++k ) {
				for( int j = 1; j <= zne.ny; ++j ) {
				for( int i = 1; i <= zne.nx; ++i )
				{
					Tcoord3D& node = zne.coordIJK(i, j, k);

					if( zne.is <= i && i <= zne.ie &&
						zne.js <= j && j <= zne.je &&
						zne.ks <= k && k <= zne.ke
						)
					{
						int n = (i - cgZne.is1) + (j - cgZne.js1)*nx0 + (k - cgZne.ks1)*nx0*ny0;
						node.x = pkXYZ[0][n];
						node.y = pkXYZ[1][n];
						node.z = pkXYZ[2][n];
					}
					else
					{
						node.x = NAN;
						node.y = NAN;
						node.z = NAN;
					}
					node.Rw = HUGE_VAL;
				}}}
			}
			else
			{
				for( int j = 1; j <= zne.ny; ++j ) {
				for( int i = 1; i <= zne.nx; ++i )
				{
					Tcoord2D& node = zne.coordIJ(i, j);

					if( zne.is <= i && i <= zne.ie &&
						zne.js <= j && j <= zne.je )
					{
						int n = (i - cgZne.is1) + (j - cgZne.js1)*nx0;
						node.x = pkXYZ[0][n];
						node.y = pkXYZ[1][n];
					}
					else
					{
						node.x = NAN;
						node.y = NAN;
					}
					node.Rw = HUGE_VAL;
				}}
			}
		}

	}  // loop through zones

	return true;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
