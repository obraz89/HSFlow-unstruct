#define _CRT_SECURE_NO_WARNINGS
// ConvertFieldUnstr2Str.cpp : Defines the entry point for the console application.
//

#include "mpi.h"

#include "optParse.h"

#include "logging.h"

#include "settings.h"

// unstruct 
#include "CGNS-ctx-unst.h"
#include "dom_uvwpt_unstruct.h"

//struct
#include "grid-load-cgns-struct.h"
#include "io-field_struct.h"

#include "interpolate_field.h"

#if defined(_WINDOWS)
  #include <windows.h>
  #include <direct.h>   // chdir
  #define chdir  _chdir

  #if ! defined(NDEBUG)
    #include <synchapi.h> // Sleep()
    #define sleep(s)   Sleep((s)*1000)
  #endif
#endif

const char szTITLE[] = "HSFlow : High-Speed Flow solver. (c) 2004-2020 HSFlow team";

TState G_State;

bool processCmdLine(int argc, char* argv[])
{
	hsLogMessage("\n\t%s\n", szTITLE);

	optparse::OptionParser parser = optparse::OptionParser()
		.usage("Usage: %prog [options] CASE_DIR")
		.version("\tBuild date: " __DATE__ " " __TIME__)
		.epilog("Positional arguments:\n"
			"  CASE_DIR\tdirectory of the case to solve");

	parser.add_option("-l", "--log").type("string").metavar("FILE").help("additionally save log to FILE");
#if ! defined(NDEBUG)
	parser.add_option("-w").type("int").metavar("N").help("wait for N seconds before start, e.g. to attach debugger");
#endif

	const optparse::Values& options = parser.parse_args(argc, argv);
	const std::vector<std::string>& args = parser.args();

	// Case dir
	if (args.size() != 1)
	{
		hsLogMessage("%s", parser.format_help().c_str());
		return false;
	}
	const std::string& dir = args.front();
	if (::chdir(dir.c_str()) < 0)
	{
		hsLogError("Can't access the case directory '%s'", dir.c_str());
		return false;
	}

	// Logging to file
	const std::string& log_fn = options["log"];
	if (!log_fn.empty())
	{
		if (!hsflow::TLog::set_file(log_fn))
			hsLogWarning("Can't open log file '%s'", log_fn.c_str());
		hsflow::TLog::log_raw_from_root(szTITLE, true);
	}

#if ! defined(NDEBUG)
	if (options.is_set("w"))
	{
		hsLogMessage("Sleeping for %d s. A debugger can be attached now...", (int)options.get("w"));
		sleep((int)options.get("w"));

#if defined(_WINDOWS)
		// Make arithmetic exceptions signalling
		unsigned int enableBits = _EM_OVERFLOW | _EM_INVALID;
		_clearfp();
		_controlfp_s(NULL, ~enableBits, enableBits);
#endif
	}
#endif

	return true;
}


int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &G_State.mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &G_State.mpiNProcs);

	hsLogMessage("My rank=%d, total procs=%d", G_State.mpiRank, G_State.mpiNProcs);

	processCmdLine(argc, argv);
	//hsLogMessage("Hello World");

	load_settings();

	// loading unstruct stuff

	G_CTXUnst.readMesh(g_Settings.strGridFnUstr);

	g_pDomUnst = &G_DomUnst;

	G_DomUnst.initializeFromCtxStage1();

	G_DomUnst.initializeFromCtxStage2();

	hsLogMessage("Calculating cell weights for vertex pv reconstruction");
	G_DomUnst.calcCellWeightsForVertices();

	G_DomUnst.allocateFlowSolution();

	G_DomUnst.initializeFlow();

	//G_DomUnst.

	// loading struct stuff

	G_Domain.nDim = 3;
	G_Domain.nu = NConsVars;

	hsLogMessage("Loading struct grid %s", g_Settings.strGridFnStrc.c_str());
	doLoadGrid_cgns_struct(g_Settings.strGridFnStrc);

	hsLogMessage("Allocating struct field");
	initField_struct();

	hsLogMessage("Interpolating field data");
	interpolate_match1to1();

	MPI_Finalize();

    return 0;
}

