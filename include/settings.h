#pragma once

#include "IniFile.hpp"

#include "logging.h"

//
// Constants
//
extern const char* g_CASE_SETTINGS_DIR;
extern const char* g_CASE_RESULTS_DIR;

// enum int with str code
// 0 <=> Str1
// 1 <=> Str2
// etc
struct t_EnumStr {

	int Val;

	std::vector<std::string> ValsStr;

	virtual void initValsStr() = 0;
	virtual const std::string& defaultValStr() const= 0;

	void set(std::string str) {
		for (int i = 0; i < ValsStr.size(); i++) {
			if (ValsStr[i].compare(str) == 0) {
				Val = i;
				return;
			}
		}
		hsLogError("EnumStr: unknown option %s", str.c_str());
	}

	std::string getOptionsStr() {
		std::string str;
		for (auto elem : ValsStr) str += elem + ",";
		return str;
	}

};

struct t_NonDimType : public t_EnumStr {
	static const int FreeStreamVelo = 0;
	static const int FreeStreamSoundSpeed = 1;
	void initValsStr() {
		ValsStr.push_back("FreeStreamVelo");
		ValsStr.push_back("FreeStreamSoundSpeed");
	}

	bool operator==(int v) { return Val == v; }

	const std::string& defaultValStr() const { return ValsStr[0]; }

	t_NonDimType() { initValsStr(); }
};

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

	// euler, ns, rans etc
	std::string strCase;

	// non-dim type
	t_NonDimType nonDimType;

	// TODO: scheme (domain) options
	double CFL;	// Courant number
	std::string strScheme;
	std::string strRiemannSolver;

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

void load_case(std::string& ini_data);

void load_case_euler(std::string& ini_data);

void load_case_ns(std::string& ini_data);