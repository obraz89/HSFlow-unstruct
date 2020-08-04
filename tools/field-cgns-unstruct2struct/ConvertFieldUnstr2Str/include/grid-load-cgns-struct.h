///////////////////////////////////////////////////////////////////////////////
// Name:        grid-load-cgns.h
// Purpose:     Definitions for CGNS grid loading
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>  // for std::lower_bound() - binary find

#include "common_data_struct.h"

typedef int PetscInt;

//-----------------------------------------------------------------------------

// CGNS family name for "frozen" zones
const char g_CG_FROZEN_ZONE_FAMILY_NAME[] = "frozen";
const int CG_MAX_NAME_LENGTH = 32 + 1/*terminating 0*/;
//-----------------------------------------------------------------------------


/**
 *  Compact storage of global indices of all nodes in a zone:
 *     Indices of ghost nodes are stored,
 *     indices of real nodes are computed
 */
class TZoneGlobIndices
{
	const TZone& zne;  // owner zone

	struct Tidx	{
		int loc;  PetscInt glob;
		bool operator<(const Tidx& that) const {  return loc < that.loc;  }
	};
	std::vector<Tidx> vec;  // NB: not using std::map to save memory

public:
	TZoneGlobIndices() = delete;
	TZoneGlobIndices(const TZone& aZone) : zne(aZone)
	{
		// Fill-in by ghost node indices
		const int nx0 = (zne.ie - zne.is + 1),
				  ny0 = (zne.je - zne.js + 1),
				  nz0 = (zne.ke - zne.ks + 1);
		vec.reserve( zne.nodes_count() - nx0*ny0*nz0  );

		for( int k=1; k<=zne.nz; ++k ){
		for( int j=1; j<=zne.ny; ++j ){
		for( int i=1; i<=zne.nx; ++i )
		{
			if( zne.is <= i && i <= zne.ie &&
				zne.js <= j && j <= zne.je &&
				zne.ks <= k && k <= zne.ke    )
				continue;

			Tidx idx = { zne.flatIdx(i,j,k), -1 };
			vec.emplace_back(idx);
		}}}
	}

	PetscInt& operator()(int i, int j, int k = 1)
	{
		// Real node: compute global index
		if( zne.is <= i && i <= zne.ie &&
			zne.js <= j && j <= zne.je &&
			zne.ks <= k && k <= zne.ke
		  )
		{
			static PetscInt ind;
			ind = zne.globRealInd(i,j,k);
			return ind;
		}

		const Tidx locIdx = { zne.flatIdx(i,j,k), -2 };
		std::vector<Tidx>::iterator it = std::lower_bound(vec.begin(), vec.end(), locIdx);
		assert( it != vec.end() && !(locIdx < *it) );
		return it->glob;
	}
};
//-----------------------------------------------------------------------------


/**
 * Temporary connectivity data for CGNS zone
 */
struct TcgnsZone
{
	enum class FaceType : char { unknown, abutted, boundary };
	FaceType face_types[6] = { };  // zero-initialization, i.e. all `unknown`

	int is1,ie1, js1,je1, ks1,ke1;  // real nodes indices range when ghosts were added, but no faces were skipped

	// Abutted faces data. The face consists of one or many patches
	struct TFacePatch
	{
		int is0,ie0, js0,je0, ks0,ke0;   // patch indices in original numbering (no ghosts, no skipped faces), may be unsored i.e. is0 > ie0
		TZoneFacePos posFace = faceNone; // face the patch belongs to

		int iZneDnr = -1;                    // 0-based index of donor zone
		TZoneFacePos dnrPosFace = faceNone;  // donor face we abutting to
		int dnr_is0, dnr_js0, dnr_ks0;  // donor patch starting indices in original numbering (no ghosts, no skipped faces)

		// Transform matrix between node indices of the current and donor zones
		int matTrans[9];
	};
	int nPatches = 0;
	TFacePatch* Patches = nullptr;

	TZoneGlobIndices* globIndices = nullptr;

	TcgnsZone() = default;
	~TcgnsZone() {  delete[] Patches;  delete globIndices;  }
};
//-----------------------------------------------------------------------------

/**
 * Info on the CGNS file to share between functions
 */
struct TcgnsContext
{
	int fileID = -1, iBase = 1;
	std::map<std::string, int>  map_ZneName2Idx;

	TcgnsZone* cgZones = nullptr;

	TcgnsContext() = default;
	~TcgnsContext(){  delete[] cgZones;  }
};
//-----------------------------------------------------------------------------
