/**********************************
// Name:        io-field.h
// Purpose:     Field input-output interfaces, cgns format
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
**********************************/

#pragma once

#include <string>

#include "common_data.h"

// Correspondence of zone number and MPI rank working with it
// g_zoneMPIRank[zone_num] -> mpi_rank
extern int* g_zoneMPIRank;

// packed zone filled by domain which is easy to handle in io operations
// TODO: make it later, now domain provides only flow
// mesh io is from cgns ctx 

static const char* g_cgCoordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };

struct t_ZonePacked {

};

//-----------------------------------------------------------------------------

bool initField();
double loadField(const std::string& fileName, const short time_layer = 0);

bool saveField(const std::string& fileName, const std::string& gridFileName,
	const short time_layer = 0);

bool read_zone_cgns(const int fileID, const int iBase, const int idxZne,
					 TpakArraysDyn<double>& field);
//-----------------------------------------------------------------------------

