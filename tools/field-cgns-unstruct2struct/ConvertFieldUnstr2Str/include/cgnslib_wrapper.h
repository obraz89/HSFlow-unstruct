///////////////////////////////////////////////////////////////////////////////
// Name:    cgnslib_wrapper.h
// Purpose: Wrapper for loading cgnslib.h with CG_ prefixes
// Author:  Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once
//-----------------------------------------------------------------------------

#include <cgnstypes.h>
#if CG_BUILD_SCOPE == 0
    #undef  CG_BUILD_SCOPE
    #define CG_BUILD_SCOPE 1
#endif

#include <cgnslib.h>
//-----------------------------------------------------------------------------
