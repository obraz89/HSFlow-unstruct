///////////////////////////////////////////////////////////////////////////////
// Project: HSFlow shared
// Purpose: Common data and functions for main program and plugins
///////////////////////////////////////////////////////////////////////////////
// File:        grid-data.cpp
// Purpose:     Implementation of TZone methods concerning 3D grid data
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#include "common_data_struct.h"
//-----------------------------------------------------------------------------


static inline const Tcoord3D& coordHalf(const Tcoord3D& n1, const Tcoord3D& n2)
{
	static Tcoord3D nH; // node half

	nH.vec = (n1.vec + n2.vec) * 0.5;
	nH.Rw = (n1.Rw + n2.Rw) * 0.5;

	return nH;
}

/**
 *  Cartesian coordinates of half-node (i+1/2, j, k)
 */
const Tcoord3D& TZone::coordI05JK(int i, int j, int k) const
{
	return ( i != nx )
		? coordHalf( coordIJK(i,j,k), coordIJK(i+1,j,k) )
		: coordIJK(i,j,k);
}

/**
 *  Cartesian coordinates of half-node (i, j+1/2, k)
 */
const Tcoord3D& TZone::coordIJ05K(int i, int j, int k) const
{
	return ( j != ny )
		? coordHalf( coordIJK(i,j,k), coordIJK(i,j+1,k) )
		: coordIJK(i,j,k);
}

/**
 *  Cartesian coordinates of half-node (i, j, k+1/2)
 */
const Tcoord3D& TZone::coordIJK05(int i, int j, int k) const
{
	return ( k != nz )
		? coordHalf( coordIJK(i,j,k), coordIJK(i,j,k+1) )
		: coordIJK(i,j,k);
}
//-----------------------------------------------------------------------------


/**
 *  Jacobian of metric coefficients
 *       / dxdk dxde dxdd \
 *  = det| dydk dyde dydd |
 *       \ dzdk dzde dzdd /
 */
inline double jacobian( const Tmtr3D& v )
{
	return
		  v.x_dksi * v.y_deta * v.z_dzet
	    - v.x_deta * v.y_dksi * v.z_dzet
	    - v.x_dksi * v.y_dzet * v.z_deta
	    + v.x_dzet * v.y_dksi * v.z_deta
	    + v.x_deta * v.y_dzet * v.z_dksi
	    - v.x_dzet * v.y_deta * v.z_dksi;
}
//-----------------------------------------------------------------------------

/**
 *  Metric coefficients in node (i,j,k) for coordinates transformation
 *  from Cartesian ("physical") to curvilinear orthogonal ("computational")
 *
 *  WARNING: Returns reference to the local static variable,
 *           contents of which change on next function call
 */
const Tmtr3D& TZone::mtrIJK(int i, int j, int k)
{
	static int prevIdx = 0;
	static Tmtr3D D;

	// Use cached data
	const int ijk = flatIdx(i,j,k);
	if( ijk == prevIdx )  // NB: assumes that on switching zones ijk also changes
		return D;

	prevIdx = ijk;
	D._invExpired = true;

	// Can't go backward
	const bool noBakI = (i == 1) || (globInd(i-1,j  ,k  ) == -1);
	const bool noBakJ = (j == 1) || (globInd(i,  j-1,k  ) == -1);
	const bool noBakK = (k == 1) || (globInd(i,  j,  k-1) == -1);

	// Can't go forward
	const bool noFwdI = (i == nx) || (globInd(i+1,j  ,k  ) == -1);
	const bool noFwdJ = (j == ny) || (globInd(i,  j+1,k  ) == -1);
	const bool noFwdK = (k == nz) || (globInd(i,  j,  k+1) == -1);

	// Grid steps in computational space
	const double l_2dksi = 0.5 * grd.l_dksi, l_2deta = 0.5 * grd.l_deta, l_2dzet = 0.5 * grd.l_dzet;

//---

	//
	// d?_dksi
	//
	if( noBakI )
	{
		const Tvec3D& v0 = coordIJK(   ijk ).vec;
		const Tvec3D& v1 = coordIJK(i+1,j,k).vec;
		const Tvec3D& v2 = coordIJK(i+2,j,k).vec;
		D.ksi = (-3.*v0 + 4.*v1 - v2) * l_2dksi;
	}
	else if( noFwdI )
	{
		const Tvec3D& v0 = coordIJK(   ijk ).vec;
		const Tvec3D& v1 = coordIJK(i-1,j,k).vec;
		const Tvec3D& v2 = coordIJK(i-2,j,k).vec;
		D.ksi = (3.*v0 - 4.*v1 + v2) * l_2dksi;
	}
	else
	{
		const Tvec3D& v1 = coordIJK(i-1,j,k).vec;
		const Tvec3D& v2 = coordIJK(i+1,j,k).vec;
		D.ksi = (v2 - v1) * l_2dksi;

		// Guard a point of tripple zones intersection
		assert( fabs(D.x_dksi) + fabs(D.y_dksi) + fabs(D.z_dksi) > 1e-10 );
	}

	//
	// d?_deta
	//
	if( noBakJ )
	{
		const Tvec3D& v0 = coordIJK(  ijk  ).vec;
		const Tvec3D& v1 = coordIJK(i,j+1,k).vec;
		const Tvec3D& v2 = coordIJK(i,j+2,k).vec;
		D.eta = (-3.*v0 + 4.*v1 - v2) * l_2deta;
	}
	else if( noFwdJ )
	{
		const Tvec3D& v0 = coordIJK(  ijk  ).vec;
		const Tvec3D& v1 = coordIJK(i,j-1,k).vec;
		const Tvec3D& v2 = coordIJK(i,j-2,k).vec;
		D.eta = (3.*v0 - 4.*v1 + v2) * l_2deta;
	}
	else
	{
		const Tvec3D& v1 = coordIJK(i,j-1,k).vec;
		const Tvec3D& v2 = coordIJK(i,j+1,k).vec;
		D.eta = (v2 - v1) * l_2deta;

		// Guard a point of tripple zones intersection
		assert( fabs(D.x_deta) + fabs(D.y_deta) + fabs(D.z_deta) > 1e-10 );
	}

	//
	// d?_dzet
	//
	if( noBakK )
	{
		const Tvec3D& v0 = coordIJK( ijk   ).vec;
		const Tvec3D& v1 = coordIJK(i,j,k+1).vec;
		const Tvec3D& v2 = coordIJK(i,j,k+2).vec;
		D.zet = (-3.*v0 + 4.*v1 - v2) * l_2dzet;
	}
	else if( noFwdK )
	{
		const Tvec3D& v0 = coordIJK( ijk   ).vec;
		const Tvec3D& v1 = coordIJK(i,j,k-1).vec;
		const Tvec3D& v2 = coordIJK(i,j,k-2).vec;
		D.zet = (3.*v0 - 4.*v1 + v2) * l_2dzet;
	}
	else
	{
		const Tvec3D& v1 = coordIJK(i,j,k-1).vec;
		const Tvec3D& v2 = coordIJK(i,j,k+1).vec;
		D.zet = (v2 - v1) * l_2dzet;

		// Guard a point of tripple zones intersection
		assert( fabs(D.x_dzet) + fabs(D.y_dzet) + fabs(D.z_dzet) > 1e-10 );
	}

	D.jac = jacobian(D);

	return D;
}
//-----------------------------------------------------------------------------


/**
 *  Metric coefficients in half-node (i+1/2,j,k) for coordinates transformation
 *  from Cartesian ("physical") to curvilinear orthogonal ("computational")
 *
 *  WARNING: Returns reference to the local static variable,
 *           contents of which change on next function call
 */
const Tmtr3D& TZone::mtrI05JK(int i, int j, int k)
{
	static int prevIdx = 0;
	static Tmtr3D D;

	// Use cached data
	const int ijk = flatIdx(i,j,k);
	if( ijk == prevIdx )  // NB: assumes that on switching zones ijk also changes
		return D;

	prevIdx = ijk;
	D._invExpired = true;

	// Can't go backward
	const bool noBakI = (i == 1) || (globInd(i-1,j  ,k  ) == -1);
	const bool noBakJ = (j == 1) || (globInd(i,  j-1,k  ) == -1);
	const bool noBakK = (k == 1) || (globInd(i,  j,  k-1) == -1);

	// Can't go forward
	const bool noFwdI = (i == nx) || (globInd(i+1,j  ,k  ) == -1);
	const bool noFwdJ = (j == ny) || (globInd(i,  j+1,k  ) == -1);
	const bool noFwdK = (k == nz) || (globInd(i,  j,  k+1) == -1);

	// Grid steps in computational space
	const double l_dksi = grd.l_dksi, l_4deta = 0.25 * grd.l_deta, l_4dzet = 0.25 * grd.l_dzet;

//---

	if( noFwdI ){
		D = mtrIJK(i,j,k);
		return D;
	}

	//
	// d?_dksi
	//
	{
		const Tvec3D& v0 = coordIJK(  ijk  ).vec;
		const Tvec3D& v1 = coordIJK(i+1,j,k).vec;
		D.ksi = (v1 - v0) * l_dksi;
	}

	//
	// d?_deta
	//
	if( noBakJ )
	{
		const Tvec3D i05j0_05 = coordIJK(   ijk   ).vec + coordIJK(i+1,j  ,k).vec;
		const Tvec3D i05j1_05 = coordIJK(i  ,j+1,k).vec + coordIJK(i+1,j+1,k).vec;
		const Tvec3D i05j2_05 = coordIJK(i  ,j+2,k).vec + coordIJK(i+1,j+2,k).vec;

		D.eta = ( -3.*i05j0_05 + 4.*i05j1_05 - i05j2_05 ) * l_4deta;
	}
	else if( noFwdJ )
	{
		const Tvec3D i05j0_05 = coordIJK(   ijk   ).vec + coordIJK(i+1,j  ,k).vec;
		const Tvec3D i05j1_05 = coordIJK(i  ,j-1,k).vec + coordIJK(i+1,j-1,k).vec;
		const Tvec3D i05j2_05 = coordIJK(i  ,j-2,k).vec + coordIJK(i+1,j-2,k).vec;

		D.eta = ( 3.*i05j0_05 - 4.*i05j1_05 + i05j2_05 ) * l_4deta;
	}
	else
	{
		const Tvec3D i05j1_05 = coordIJK(i  ,j-1,k).vec + coordIJK(i+1,j-1,k).vec;
		const Tvec3D i05j2_05 = coordIJK(i  ,j+1,k).vec + coordIJK(i+1,j+1,k).vec;

		D.eta = ( i05j2_05 - i05j1_05 ) * l_4deta;
	}

	//
	// d?_dzet
	//
	if( noBakK )
	{
		const Tvec3D i05k0_05 = coordIJK(  ijk    ).vec + coordIJK(i+1,j,k  ).vec;
		const Tvec3D i05k1_05 = coordIJK(i  ,j,k+1).vec + coordIJK(i+1,j,k+1).vec;
		const Tvec3D i05k2_05 = coordIJK(i  ,j,k+2).vec + coordIJK(i+1,j,k+2).vec;

		D.zet = ( -3.*i05k0_05 + 4.*i05k1_05 - i05k2_05 ) * l_4dzet;
	}
	else if( noFwdK )
	{
		const Tvec3D i05k0_05 = coordIJK(   ijk   ).vec + coordIJK(i+1,j,k  ).vec;
		const Tvec3D i05k1_05 = coordIJK(i  ,j,k-1).vec + coordIJK(i+1,j,k-1).vec;
		const Tvec3D i05k2_05 = coordIJK(i  ,j,k-2).vec + coordIJK(i+1,j,k-2).vec;

		D.zet = ( 3.*i05k0_05 - 4.*i05k1_05 + i05k2_05 ) * l_4dzet;
	}
	else
	{
		const Tvec3D i05k1_05 = coordIJK(i  ,j,k-1).vec + coordIJK(i+1,j,k-1).vec;
		const Tvec3D i05k2_05 = coordIJK(i  ,j,k+1).vec + coordIJK(i+1,j,k+1).vec;

		D.zet = ( i05k2_05 - i05k1_05 ) * l_4dzet;
	}

	D.jac = jacobian(D);

	return D;
}
//-----------------------------------------------------------------------------


/**
 *  Metric coefficients in half-node (i,j+1/2,k) for coordinates transformation
 *  from Cartesian ("physical") to curvilinear orthogonal ("computational")
 *
 *  WARNING: Returns reference to the local static variable,
 *           contents of which change on next function call
 */
const Tmtr3D& TZone::mtrIJ05K(int i, int j, int k)
{
	static int prevIdx = 0;
	static Tmtr3D D;

	// Use cached data
	const int ijk = flatIdx(i,j,k);
	if( ijk == prevIdx )  // NB: assumes that on switching zones ijk also changes
		return D;

	prevIdx = ijk;
	D._invExpired = true;

	// Can't go backward
	const bool noBakI = (i == 1) || (globInd(i-1,j  ,k  ) == -1);
	const bool noBakJ = (j == 1) || (globInd(i,  j-1,k  ) == -1);
	const bool noBakK = (k == 1) || (globInd(i,  j,  k-1) == -1);

	// Can't go forward
	const bool noFwdI = (i == nx) || (globInd(i+1,j  ,k  ) == -1);
	const bool noFwdJ = (j == ny) || (globInd(i,  j+1,k  ) == -1);
	const bool noFwdK = (k == nz) || (globInd(i,  j,  k+1) == -1);

	// Grid steps in computational space
	const double l_4dksi = 0.25 * grd.l_dksi, l_deta = grd.l_deta, l_4dzet = 0.25 * grd.l_dzet;

//---

	if( noFwdJ ){
		D = mtrIJK(i,j,k);
		return D;
	}

	//
	// d?_dksi
	//
	if( noBakI )
	{
		const Tvec3D j05i0_05 = coordIJK(   ijk   ).vec + coordIJK(i  ,j+1,k).vec;
		const Tvec3D j05i1_05 = coordIJK(i+1,j  ,k).vec + coordIJK(i+1,j+1,k).vec;
		const Tvec3D j05i2_05 = coordIJK(i+2,j  ,k).vec + coordIJK(i+2,j+1,k).vec;

		D.ksi = ( -3.*j05i0_05 + 4.*j05i1_05 - j05i2_05 ) * l_4dksi;
	}
	else if( noFwdI )
	{
		const Tvec3D j05i0_05 = coordIJK(   ijk   ).vec + coordIJK(i  ,j+1,k).vec;
		const Tvec3D j05i1_05 = coordIJK(i-1,j  ,k).vec + coordIJK(i-1,j+1,k).vec;
		const Tvec3D j05i2_05 = coordIJK(i-2,j  ,k).vec + coordIJK(i-2,j+1,k).vec;

		D.ksi = ( 3.*j05i0_05 - 4.*j05i1_05 + j05i2_05 ) * l_4dksi;
	}
	else
	{
		const Tvec3D j05i1_05 = coordIJK(i-1,j  ,k).vec + coordIJK(i-1,j+1,k).vec;
		const Tvec3D j05i2_05 = coordIJK(i+1,j  ,k).vec + coordIJK(i+1,j+1,k).vec;

		D.ksi = ( j05i2_05 - j05i1_05 ) * l_4dksi;
	}

	//
	// d?_deta
	//
	{
		const Tvec3D& v0 = coordIJK(i  ,j  ,k).vec;
		const Tvec3D& v1 = coordIJK(i  ,j+1,k).vec;

		D.eta = ( v1 - v0 ) * l_deta;
	}

	//
	// d?_dzet
	//
	if( noBakK )
	{
		const Tvec3D j05k0_05 = coordIJK( ijk     ).vec + coordIJK(i,j+1,k  ).vec;
		const Tvec3D j05k1_05 = coordIJK(i,j  ,k+1).vec + coordIJK(i,j+1,k+1).vec;
		const Tvec3D j05k2_05 = coordIJK(i,j  ,k+2).vec + coordIJK(i,j+1,k+2).vec;

		D.zet = ( -3.*j05k0_05 + 4.*j05k1_05 - j05k2_05 ) * l_4dzet;
	}
	else if( noFwdK )
	{
		const Tvec3D j05k0_05 = coordIJK( ijk     ).vec + coordIJK(i,j+1,k  ).vec;
		const Tvec3D j05k1_05 = coordIJK(i,j  ,k-1).vec + coordIJK(i,j+1,k-1).vec;
		const Tvec3D j05k2_05 = coordIJK(i,j  ,k-2).vec + coordIJK(i,j+1,k-2).vec;

		D.zet = ( 3.*j05k0_05 - 4.*j05k1_05 + j05k2_05 ) * l_4dzet;
	}
	else
	{
		const Tvec3D j05k1_05 = coordIJK(i,j  ,k-1).vec + coordIJK(i,j+1,k-1).vec;
		const Tvec3D j05k2_05 = coordIJK(i,j  ,k+1).vec + coordIJK(i,j+1,k+1).vec;

		D.zet = ( j05k2_05 - j05k1_05 ) * l_4dzet;
	}

	D.jac = jacobian(D);

	return D;
}
//-----------------------------------------------------------------------------


/**
 *  Metric coefficients in half-node (i,j,k+1/2) for coordinates transformation
 *  from Cartesian ("physical") to curvilinear orthogonal ("computational")
 *
 *  WARNING: Returns reference to the local static variable,
 *           contents of which change on next function call
 */
const Tmtr3D& TZone::mtrIJK05(int i, int j, int k)
{
	static int prevIdx = 0;
	static Tmtr3D D;

	// Use cached data
	const int ijk = flatIdx(i,j,k);
	if( ijk == prevIdx )  // NB: assumes that on switching zones ijk also changes
		return D;

	prevIdx = ijk;
	D._invExpired = true;

	// Can't go backward
	const bool noBakI = (i == 1) || (globInd(i-1,j  ,k  ) == -1);
	const bool noBakJ = (j == 1) || (globInd(i,  j-1,k  ) == -1);
	const bool noBakK = (k == 1) || (globInd(i,  j,  k-1) == -1);

	// Can't go forward
	const bool noFwdI = (i == nx) || (globInd(i+1,j  ,k  ) == -1);
	const bool noFwdJ = (j == ny) || (globInd(i,  j+1,k  ) == -1);
	const bool noFwdK = (k == nz) || (globInd(i,  j,  k+1) == -1);

	// Grid steps in computational space
	const double l_4dksi = 0.25 * grd.l_dksi, l_4deta = 0.25 * grd.l_deta, l_dzet = grd.l_dzet;

//---

	if( noFwdK ){
		D = mtrIJK(i,j,k);
		return D;
	}

	//
	// d?_dksi
	//
	if( noBakI )
	{
		const Tvec3D k05i0_05 = coordIJK(   ijk   ).vec + coordIJK(i  ,j,k+1).vec;
		const Tvec3D k05i1_05 = coordIJK(i+1,j,k  ).vec + coordIJK(i+1,j,k+1).vec;
		const Tvec3D k05i2_05 = coordIJK(i+2,j,k  ).vec + coordIJK(i+2,j,k+1).vec;

		D.ksi = ( -3.*k05i0_05 + 4.*k05i1_05 - k05i2_05 ) * l_4dksi;
	}
	else if( noFwdI )
	{
		const Tvec3D k05i0_05 = coordIJK(   ijk   ).vec + coordIJK(i  ,j,k+1).vec;
		const Tvec3D k05i1_05 = coordIJK(i-1,j,k  ).vec + coordIJK(i-1,j,k+1).vec;
		const Tvec3D k05i2_05 = coordIJK(i-2,j,k  ).vec + coordIJK(i-2,j,k+1).vec;

		D.ksi = ( 3.*k05i0_05 - 4.*k05i1_05 + k05i2_05 ) * l_4dksi;
	}
	else
	{
		const Tvec3D k05i1_05 = coordIJK(i-1,j,k  ).vec + coordIJK(i-1,j,k+1).vec;
		const Tvec3D k05i2_05 = coordIJK(i+1,j,k  ).vec + coordIJK(i+1,j,k+1).vec;

		D.ksi = ( k05i2_05 - k05i1_05 ) * l_4dksi;
	}

	//
	// d?_deta
	//
	if( noBakJ )
	{
		const Tvec3D k05j0_05 = coordIJK( ijk     ).vec + coordIJK(i,j  ,k+1).vec;
		const Tvec3D k05j1_05 = coordIJK(i,j+1,k  ).vec + coordIJK(i,j+1,k+1).vec;
		const Tvec3D k05j2_05 = coordIJK(i,j+2,k  ).vec + coordIJK(i,j+2,k+1).vec;

		D.eta = ( -3.*k05j0_05 + 4.*k05j1_05 - k05j2_05 ) * l_4deta;
	}
	else if( noFwdJ )
	{
		const Tvec3D k05j0_05 = coordIJK( ijk     ).vec + coordIJK(i,j  ,k+1).vec;
		const Tvec3D k05j1_05 = coordIJK(i,j-1,k  ).vec + coordIJK(i,j-1,k+1).vec;
		const Tvec3D k05j2_05 = coordIJK(i,j-2,k  ).vec + coordIJK(i,j-2,k+1).vec;

		D.eta = ( 3.*k05j0_05 - 4.*k05j1_05 + k05j2_05 ) * l_4deta;
	}
	else
	{
		const Tvec3D k05j1_05 = coordIJK(i,j-1,k  ).vec + coordIJK(i,j-1,k+1).vec;
		const Tvec3D k05j2_05 = coordIJK(i,j+1,k  ).vec + coordIJK(i,j+1,k+1).vec;

		D.eta = ( k05j2_05 - k05j1_05 ) * l_4deta;
	}

	//
	// d?_dzet
	//
	{
		const Tvec3D& v0 = coordIJK( ijk   ).vec;
		const Tvec3D& v1 = coordIJK(i,j,k+1).vec;
		D.zet = (v1 - v0) * l_dzet;
	}

	D.jac = jacobian(D);

	return D;
}
//-----------------------------------------------------------------------------


const Tmtr3D::inv& Tmtr3D::getInverseMetric() const
{
	if( _invExpired )
	{
		// Direct and Inverse metric coefficients transformation:
		//  | dxdk dxde dxdd |   | dkdx dkdy dkdz |   | 1 0 0 |
		//  | dydk dyde dydd | x | dedx dedy dedz | = | 0 1 0 |
		//  | dzdk dzde dzdd |   | dddx dddy dddz |   | 0 0 1 |
		assert( fabs(jac) > 1e-16 );
		inv& D = _inv;

		D.ksi_dx = ( y_deta*z_dzet - y_dzet*z_deta)/jac;
		D.ksi_dy = (-x_deta*z_dzet + x_dzet*z_deta)/jac;
		D.ksi_dz = ( x_deta*y_dzet - x_dzet*y_deta)/jac;

		D.eta_dx = (-y_dksi*z_dzet + y_dzet*z_dksi)/jac;
		D.eta_dy = ( x_dksi*z_dzet - x_dzet*z_dksi)/jac;
		D.eta_dz = (-x_dksi*y_dzet + x_dzet*y_dksi)/jac;

		D.zet_dx = ( y_dksi*z_deta - y_deta*z_dksi)/jac;
		D.zet_dy = (-x_dksi*z_deta + x_deta*z_dksi)/jac;
		D.zet_dz = ( x_dksi*y_deta - x_deta*y_dksi)/jac;

		_invExpired = false;
	}

	return _inv;
}
//-----------------------------------------------------------------------------
