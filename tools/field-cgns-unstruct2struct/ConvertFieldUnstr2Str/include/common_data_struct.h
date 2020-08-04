///////////////////////////////////////////////////////////////////////////////
// Name:        common_data.h
// Purpose:     Common data shared through main exe and plugins
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include <map>
#include <assert.h>
#include <type_traits>

#include "common_data.h"

#include "evector.h"

#pragma warning(push)
#pragma warning(disable:4251)  // class 'std::string' needs to have dll-interface

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

//
// Grid coordinates in Cartesian ("physical") space
//
template <bool is3D> struct Tcoord;

template<>
struct Tcoord<false>
{
	union {
		struct { double x, y; };
		Tvec<false> vec;
	};
	double Rw;   // distance to the nearest wall
};

template<>
struct Tcoord<true>
{
	union {
		struct { double x, y, z; };
		Tvec<true> vec;
	};
	double Rw;   // distance to the nearest wall
};

using Tcoord3D = Tcoord<true>;
using Tcoord2D = Tcoord<false>;
//-----------------------------------------------------------------------------

//
// Metric coefficients of coordinates transformation
//
// (x,y,z)       - Cartesian ("physical") coordinates
// (ksi,eta,zeta) - curvilinear ("computational") coordinates
//
template <bool is3D> struct Tmtr;

template <>
struct Tmtr<false>
{
	// Direct 2D metric coefficients d{x,y}/d{ksi,eta}
	union { struct { double x_dksi, y_dksi; }; Tvec2D ksi; };
	union { struct { double x_deta, y_deta; }; Tvec2D eta; };

	double jac;  // transformation Jacobian det[d(x,y)/d(ksi,eta)]

	// Inverse 2D metric coefficients d{ksi,eta}/d{x,y}
	struct inv
	{
		union { struct { double ksi_dx, ksi_dy; }; Tvec2D ksi; };
		union { struct { double eta_dx, eta_dy; }; Tvec2D eta; };
	};
	const inv& getInverseMetric() const;

private:
	mutable inv  _inv;         // cached inverse metric data
	mutable bool _invExpired;  // cache is expired and should be recalculated

	friend struct TZone;
};
//-----------------------------------------------------------------------------

template <>
struct Tmtr<true>
{
	// Direct 3D metric coefficients d{x,y,z}/d{ksi,eta,zeta}
	union { struct { double x_dksi, y_dksi, z_dksi; }; Tvec3D ksi; };
	union { struct { double x_deta, y_deta, z_deta; }; Tvec3D eta; };
	union { struct { double x_dzet, y_dzet, z_dzet; }; Tvec3D zet; };

	double jac;  // transformation Jacobian

	// Inverse 3D metric coefficients d{ksi,eta}/d{x,y}
	struct inv
	{
		union { struct { double ksi_dx, ksi_dy, ksi_dz; }; Tvec3D ksi; };
		union { struct { double eta_dx, eta_dy, eta_dz; }; Tvec3D eta; };
		union { struct { double zet_dx, zet_dy, zet_dz; }; Tvec3D zet; };
	};
	const inv& getInverseMetric() const;

private:
	mutable inv  _inv;         // cached inverse metric data
	mutable bool _invExpired;  // cache is expired and should be recalculated

	friend struct TZone;
};

using Tmtr3D = Tmtr<true>;
using Tmtr2D = Tmtr<false>;
//-----------------------------------------------------------------------------


//
// Domain zones (aka blocks)
//
struct TZone;

struct TZoneFace
{
	// Boundary condition on the face
	char szBC[33] = "";                 // BC-family name, MUST be empty if abutted

	bool isSkipped = false;     // face's grid layer skipped for processing by abutted zone

	// Transformation for periodic BC to apply at ghost nodes
	struct Ttransform
	{
		int is,ie, js,je, ks,ke;  // indices range of nodes to transform (ghosts)

		// Transform from the donor face to the current, i.e. reversed (!!!)
		// Parallel translation vector
		float Translation[3];
		// Affine transform: rotation or reflection
		float TransformCenter[3];  double TransformMatrix[3*3];
	};
	Ttransform* transform = nullptr;
};

// Zone face position in curvilinear coords
enum TZoneFacePos : char { faceNone = -1, faceXmin = 0, faceXmax, faceYmin, faceYmax, faceZmin, faceZmax };

struct TZone
{
	char szName[40]{ '\0' };  // name of the zone, initialized by '\0'
	bool isFrozen = false; // don't compute anything in this zone, keep field fixed

	// Zone dimensions including ghost nodes
	int nx = 0, ny = 0, nz = 1;
	int nodes_count() const {
		return nx*ny*nz;
	}

//
// Indices
//
	// Start & end 1-based indices of real nodes
	int is=0,ie=0,  js=0,je=0,  ks=0,ke=0;

	// Inter-zone connectivity data
	int nGlobStart = 0;         // starting global 0-based index of the real nodes
	int* globIndices = nullptr; // global 0-based indices of each zone's node (real & ghost)
	// NB: maximum nodes count is 2'147 mln !!!

	// Local to the zone 1-based flattened 1D index of (i,j,k)-node
	int flatIdx(int i, int j) const
	{
		assert( nz <= 1 );
		//assert( i>=1 && i<=nx && j>=1 && j<=ny );
		return i + (j-1)*nx;
	}
	int flatIdx(int i, int j, int k) const
	{
		assert( nz > 1 || k == 1 );
		//assert( i>=1 && i<=nx && j>=1 && j<=ny && k>=1 && k<=nz );
		return i + (j-1)*nx + (k-1)*nx*ny;
	}

	// Global (inter-zones) 0-based index of the *real* (i,j,k)-node
	int globRealInd(int i, int j) const
	{
		assert( nz<=1 && i>=is && i<=ie && j>=js && j<=je );
		// size excluding ghosts
		const int nx0 = ie - is + 1;
		return  nGlobStart + (i-is) + (j-js)*nx0;
	}
	int globRealInd(int i, int j, int k) const
	{
		assert( i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke );
		// size excluding ghosts
		const int nx0 = ie - is + 1,  ny0 = je - js + 1;
		return  nGlobStart + (i-is) + (j-js)*nx0 + (k-ks)*nx0*ny0;
	}

	// Global (inter-zones) 0-based index of *any* (real or ghost) (i,j,k)-node
	int& globInd(int idx)
	{
		assert( globIndices && idx >= 1 && idx <= nx*ny*nz );
		return globIndices[idx - 1];
	}
	int& globInd(int i, int j)       {   return globInd( flatIdx(i,j)   );   }
	int& globInd(int i, int j, int k){   return globInd( flatIdx(i,j,k) );   }

//
// Grid data
//
	struct {
		Tcoord2D* c2d = nullptr;
		Tcoord3D* c3d = nullptr;

		// grid steps in computational space
		double dksi, deta, dzet;
		double l_dksi, l_deta, l_dzet;
	} grd;

	//
	// Coordinates in Cartesian space
	//
	Tcoord2D& coordIJ(int idx /**< 1-based flattened node index */) const {
		assert( idx > 0 && idx <= flatIdx(nx,ny) );
		return grd.c2d[ idx - 1 ];
	}
	Tcoord2D& coordIJ(int i, int j) const {         // (i,     j    )
		return coordIJ( flatIdx(i,j) );
	}
	const Tcoord2D& coordI05J(int i, int j) const;  // (i+1/2, j    )
	const Tcoord2D& coordIJ05(int i, int j) const;  // (i,     j+1/2)


	Tcoord3D& coordIJK(int idx /**< 1-based flattened node index */) const {
		assert( idx > 0 && idx <= flatIdx(nx,ny,nz) );
		return grd.c3d[ idx - 1 ];
	}
	Tcoord3D& coordIJK(int i, int j, int k) const {         // (i,     j,     k    )
		return coordIJK( flatIdx(i,j,k) );
	}
	const Tcoord3D& coordI05JK(int i, int j, int k) const;  // (i+1/2, j,     k    )
	const Tcoord3D& coordIJ05K(int i, int j, int k) const;  // (i,     j+1/2, k    )
	const Tcoord3D& coordIJK05(int i, int j, int k) const;  // (i,     j,     k+1/2)

	// Overridden methods 
	template<bool is3D, class = typename std::enable_if<is3D>::type >
	const Tcoord3D& coord(int i, int j, int k) {
		return coordIJK(i, j, k);
	}

	template<bool is3D, class = typename std::enable_if<!is3D>::type >
	const Tcoord2D& coord(int i, int j, int k) {
		assert(k <= 1);
		return coordIJ(i, j);
	}


	//
	// Metric coefficients
	//
	const Tmtr2D& mtrIJ(int i, int j);    // (i,     j    )
	const Tmtr2D& mtrI05J(int i, int j);  // (i+1/2, j    )
	const Tmtr2D& mtrIJ05(int i, int j);  // (i,     j+1/2)

	const Tmtr3D& mtrIJK(int i, int j, int k);    // (i,     j,     k    )
	const Tmtr3D& mtrI05JK(int i, int j, int k);  // (i+1/2, j,     k    )
	const Tmtr3D& mtrIJ05K(int i, int j, int k);  // (i,     j+1/2, k    )
	const Tmtr3D& mtrIJK05(int i, int j, int k);  // (i,     j,     k+1/2)

	template<bool is3D, class = typename std::enable_if<is3D>::type >
	const Tmtr3D& mtr(int i, int j, int k) {
		return mtrIJK(i, j, k);
	}

	template<bool is3D, class = typename std::enable_if<!is3D>::type >
	const Tmtr2D& mtr(int i, int j, int k) {
		return mtrIJ(i, j);
	}

	/**
	 * Linear scale of a grid cell
	 *
	 * @param mtr - metric coefficients of the cell
	 *
	 * @retval square root of the approximate area of the 2D cell
	 * @retval cubic root of the approximate volume of the 3D cell
	 */
	template<bool is3D>
	double cellScale(const Tmtr<is3D>& mtr) const
	{
		const double h = 1e-8;  // minimum scale
		if( ! is3D ) {
			const double s = mtr.jac * (grd.dksi * grd.deta);
			if( s > h*h )
				return sqrt(s);
		}
		else {
			const double v = mtr.jac * (grd.dksi * grd.deta * grd.dzet);
			if( v > h*h*h )
				return cbrt(v);
		}

		return h;
	}

	//
	// Primitive flow variables in all nodes of the zone (including ghosts)
	//
	union {
		struct {
			double* __restrict U;    // current (n+1) time layer
			double* __restrict Un;   // previous (n) time layer
			double* __restrict Unm1; // pre-previous (n-1) time layer
		};
		double* UU[3];
	};

	// Zone faces info
	TZoneFace Faces[6];  // Imin,Imax, Jmin,Jmax, Kmin,Kmax

	TZone() : U(nullptr), Un(nullptr), Unm1(nullptr) {
		;
	}

	TZone(const TZone&) = delete;
	void operator=(const TZone&) = delete;

	~TZone()
	{
		delete[] globIndices;

		delete[] U;
		delete[] Un;
		delete[] Unm1;

		delete[] grd.c2d;
		delete[] grd.c3d;
	}
};
//-----------------------------------------------------------------------------

//
// The whole computational domain
//
struct TDomain
{
	int nu,	  // number of dependent (unknown) variables
	    nDim; // number of independent variables (problem dimensions)

	// Domain zones with grid and solution data
	//
	int nZones;  // total number of zones
	int bs, be;  // start & end (inclusive) 0-based zone (aka block) indices in the current MPI rank
	int* map_iZne2cgID;  // map_iZne2cgID[b] == cgZne, where b -- internal 0-based zone index, cgZne -- CGNS 1-based zone ID
	TZone* Zones;

	// Global grid info
	double gridCellScaleMin, gridCellScaleMax;

	// Gas parameters
	double (*pfunViscosity)(const double&) = nullptr;

	// Info for input-output
	std::map<std::string, double> mapCasePrms_real;
};

bool assignZonesToProcs();

extern TDomain G_Domain;
//-----------------------------------------------------------------------------

#pragma warning(pop)
