#define _CRT_SECURE_NO_WARNINGS

#include "CGNS-ctx.h"

#include "logging.h"

#include "optParse.h"

#include "common_data.h"

#include "settings.h"

#include "mesh.h"

#include "bc_common.h"

// model-specific part
#include "bc_euler.h"
#include "flow_model_perfect_gas.h"
#include "dom_euler.h"

#if defined(_WINDOWS)
	#include <direct.h>   // chdir
	#define chdir  _chdir

	#include <winsock.h>  // gethostname

	#if ! defined(NDEBUG)
		#include <synchapi.h> // Sleep()
		#define sleep(s)   Sleep((s)*1000)
	#endif
#else
	#include <unistd.h>  // sleep
#endif

const char szTITLE[] = "HSFlow: High-Speed Flow solver. (c) 2004-2020 HSFlow team";

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
	if (G_State.mpiNProcs > 1 && options.is_set("w"))
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

int main(int argc, char* argv[])
{

	int err = EXIT_FAILURE;

	if (!processCmdLine(argc, argv))
		goto fin;

	// Solver run identification: current date & time + hostname
	{
		time_t t = time(nullptr);
		char tm[32];   strftime(tm, 32, "%Y-%m-%d %H:%M:%S", localtime(&t));

		char host[256] = {};  gethostname(host, 256);

		hsLogMessage(
			"\t-\n"
			"\t- Run %s @ %s\n"
			"\t-\n",
			tm, host
		);
	}
	//hsLogWTime(true);

	G_pMesh = &G_Domain;
	G_pBCList = &G_BCListEuler;

	if (!load_settings())
		goto fin;

	G_CGNSCtx.readMesh(g_genOpts.strGridFN);

	G_Domain.initializeFromCtx();

	G_Domain.allocateFlowSolution();

	initialize_flow_model();

	G_Domain.initializeFlow();

	G_Domain.dump_flow();

	G_Domain.makeTimeStep();

	err = EXIT_SUCCESS;

fin:
	//hsflow::TLog::destroy();  // flushes remaining messages
	delete[] G_Domain.Zones;
	return err;
}
