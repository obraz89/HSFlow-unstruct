///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from EXE
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#pragma once

// SHARED_EXPORT should be defined only in HSFlow.exe

#if defined(_WIN32)
#	if defined(SHARED_EXPORT)
#		define DLLIMPEXP  __declspec(dllexport)
#	else
#		define DLLIMPEXP  __declspec(dllimport)
#	endif
#else // linux
#	define DLLIMPEXP
#endif

#if defined(_WIN32)
#	define DLLEXPORT  __declspec(dllexport)
#else // linux
#	define DLLEXPORT
#endif

//-----------------------------------------------------------------------------
