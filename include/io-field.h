/**********************************
// Name:        io-field.h
// Purpose:     Field input-output interfaces, cgns format
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
**********************************/

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

bool read_zone_cgns(const int fileID, const int iBase, const int idxZne,
					 TpakArraysDyn<double>& grid, TpakArraysDyn<double>& field);
//-----------------------------------------------------------------------------

