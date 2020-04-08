///////////////////////////////////////////////////////////////////////////////
// Name:        common_procs.h
// Purpose:     Generic procedures shared through main exe and plugins
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <ctime>
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
// File system utilities
//
bool hs_file_exists(const std::string& fn);
time_t hs_file_mtime(const std::string& fn);
bool hs_dir_exists(const std::string& dn);
bool hs_dir_create(const std::string& dn);
//-----------------------------------------------------------------------------

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
