///////////////////////////////////////////////////////////////////////////////
// Name:        common_procs.h
// Purpose:     Generic procedures shared through main exe and plugins
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <ctime>

#include "dll_import-export.h"
//-----------------------------------------------------------------------------

//
// HSFlow logging
//
#include "logging.h"
//-----------------------------------------------------------------------------

#include "common_data.h"

/// Reimplementation std::copy_n() for older compilers not supporting it
template<typename T>
inline void copy_n(const T* from, size_t count, T* to)
{
	std::copy(from, from+count, to);
}

//
// Math functions
//

/// Squared value
template<typename T>
inline T sq(T x)
{
	return x*x;
}

// geometry routines
void ComputeTriangleAreaNormal(const t_Vec3 (&pnts)[3], t_Vec3& norm, double& area);

void ComputeQuadAreaNormal(const t_Vec3(&pnts)[4], t_Vec3& norm, double& area);
//-----------------------------------------------------------------------------
