///////////////////////////////////////////////////////////////////////////////
// Name:        io.h
// Purpose:     Input/Output functions
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

// Correspondence of zone number and MPI rank working with it
// g_zoneMPIRank[zone_num] -> mpi_rank
extern int* g_zoneMPIRank;
//-----------------------------------------------------------------------------

bool initField();
double loadField(const std::string& fileName, const short time_layer = 0);

bool saveField(const std::string& fileName, const std::string& gridFileName,
			   const short time_layer = 0, bool isDouble = true);
bool saveField_legacy(const char* aFileName);

bool saveFieldWall(const std::string& fileName, const std::string& gridFileName);
bool saveFieldSlices(const std::string& fileName_base, const std::string& gridFileName);
//-----------------------------------------------------------------------------
