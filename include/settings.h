#pragma once

#include "IniFile.hpp"

//
// Constants
//
extern const char* g_CASE_SETTINGS_DIR;
extern const char* g_CASE_RESULTS_DIR;

struct TgenericSettings
{
	//--- Init ---
	std::string strInitFieldFN; // init field filename
//	std::string strPrevFieldFN; // field at previous time step, approximating if empty
	std::string strGridFN;      // grid file name
//	bool allowFrozenZones;      // enable frozen zones if specified in the grid
//	float init_ltzero_filter; // filter out non-physical negative values replacing by `init_ltzero_filter * U_inf`, <=0 - if not needed

	//--- Time ---
	double timeStart;  // initial time, if <0 - take from init-file
	int numTimeSteps; // number of time steps to do
	double convergenceResidual; // stop solver at this residual

	double CFL;	// Courant number

	// Automatic increasing of time step size:
//	struct {
//		double factor;       // dt *= factor
//		int nSteps_interval; // check after this number of steps
//		double limit;        // increase timestep up to this value
//		double minResidual;  // increase timestep only if 0-th NES residual is lower
//	} dt_inc;

	//--- Iterations ---
//	double epsOut;   //Newton iteration residual
//	int iterOut;     //max num of Newton iterations
//	double epsIn;    //internal iterations residual (GMRES)
//	int iterIn;      //max num of internal iterations (GMRES)
//	int itsFalseJac; // max num of iterations with False Jacoby matrix
//	double startNwtM;   //init value for tau param in modified Newton method
//	double qmax;		//relative decrement of residual at Newton step -- indication to recal Jac matrix

	//--- Saving ---
	std::string strSaveFieldFNTemplate;  // printf-style template for generating file name using time
	std::string strSaveGridFN;           // name of separate grid file, if empty grid is embedded

	int timeSteps2Write;            // time steps count to save full field
	//bool isDoublePrecision2Write;   // precision to write full field

//	int slices_timeSteps2Write;    // time steps count to write field slices (<=0 if not saving)
//	std::vector<TSliceZneDef> lstSlices;
	//bool slices_isDoublePrecision2Write;  // precision to write slices

//	int wall_timeSteps2Write;          // time steps count to write wall data (<=0 if not saving)
//	std::set<std::string> setWallFuncs;

//  bool wall_isDoublePrecision2Write;   // precision to write wall
};

extern TgenericSettings g_genOpts;
//-----------------------------------------------------------------------------

bool load_settings();
