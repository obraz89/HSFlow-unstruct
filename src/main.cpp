#define _CRT_SECURE_NO_WARNINGS

#include "CGNS-ctx.h"

#include "logging.h"

#include "optParse.h"

#include "common_data.h"

#include "settings.h"

#include "dom_base.h"

#include "bc_common.h"

#include "ghost_common.h"

#include "io-field.h"

#include "petsc.h"

#include "flow_model_perfect_gas.h"

#include <ctime>

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

int main(int argc, char* argv[])
{

	int err = EXIT_FAILURE;

	//if (MPI_Init(&argc, &argv) != 0)
	if (PetscInitializeNoArguments() != 0)
	{
		printf("Can't initialize Petsc" "\n");  // don't use `hsLog*()`
		return err;
	}

	//PetscPopSignalHandler();  // don't let PETSc handle signals

	MPI_Comm_rank(MPI_COMM_WORLD, &G_State.mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &G_State.mpiNProcs);

	if (!processCmdLine(argc, argv))
		goto fin;

	hsLogMsgAllRanks("My rank=%d, total procs=%d", G_State.mpiRank, G_State.mpiNProcs);

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
	hsLogWTime(true);
	
	if (!load_settings())
		goto fin;

	G_CGNSCtx.readMesh(g_genOpts.strGridFN);

	// load cells & vertices
	// make vertex connectivity so ghost manager is able to do some funcs
	G_pDom->initializeFromCtxStage1();

	// initialize ghost manager, it uses some basic funcs of mesh
	// requires first stage of mesh init
	G_pGhostMngBase->initialize(G_CGNSCtx);

	// cell connectivity, face list
	G_pDom->initializeFromCtxStage2();


	MPI_Barrier(MPI_COMM_WORLD);
	hsflow::TLog::flush();

	// update cell center etc for ghosts
	// cuurently not required as we store full mesh at each worker 
	//G_GhostMngEu.exchangeGeomData();

	G_pDom->allocateFlowSolution();

	initialize_flow_model();

	G_pDom->initializeFlow();

	G_pDom->checkMesh();

	hsLogMessage("Computing reconstruction data...");

	G_pDom->prepareBeforeTimeMarch();

	hsflow::TLog::flush();
	{
		int count = 0; 

		for (int iTStep = 0; iTStep < g_genOpts.numTimeSteps; iTStep++) {

			hsLogMessage("Iter #%d:", iTStep);

			G_pDom->makeTimeStep();
			G_pDom->checkMinMaxCSV();

			if (++count >= g_genOpts.timeSteps2Write) {
				// TODO: G_Domain.saveFiled()
				char fld_name[64];
				sprintf(fld_name, "%.5f", G_State.time);
				bool ok = saveField(fld_name, "");
				if (!ok) hsLogError("There was an error in writing field to output file");
				count = 0;
			}

			hsflow::TLog::flush();

		}
	}

	err = EXIT_SUCCESS;

fin:
	// DEBUG: testing offsets 
	//std::vector<t_CellKindRange> offsets = G_Domain.Zones[0].getCellsOffsets();

	hsflow::TLog::destroy();  // flushes remaining messages

	//MPI_Finalize();
	PetscFinalize();

	return err;


}
