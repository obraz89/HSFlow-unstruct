#pragma once

//
// Problem solving state
//
struct TState
{
	int mpiRank, mpiNProcs;  // MPI rank, number of procs
							 // Relation of zone index and MPI rank working with it
	int* map_zone2rank;  // map_zone2rank[izne] == mpiRank, where izne -- 0-based zone index

						 //---

						 //int nTmStep;     // current time step number
						 //int nwtIter;     // current Newton's iteration number

						 /// L_inf norm of residual at 0-th Newton's iteration at current step
						 //double initialResidual;

						 // TODO: these are scheme parameters, separate them later

	double time;

	double ResidTot;

	// debug vars
	double ResidNormVeloWall;
};

bool hs_file_exists(const std::string& fn);
//time_t hs_file_mtime(const std::string& fn);
bool hs_dir_exists(const std::string& dn);
bool hs_dir_create(const std::string& dn);

// initialized in main.cpp
extern TState G_State;
