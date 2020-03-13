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
//-----------------------------------------------------------------------------
