#include "settings.h"

#include "flow_params.h"

#include "flow_model_perfect_gas.h"

void load_case_ns(std::string& ini_data) {

	G_FreeStreamParams.init(ini_data, "");

	//G_pBCList = &G_BCListEuler;

	//G_BCListEuler.init(ini_data, "");

	//G_pGhostMngBase = &G_GhostMngEu;
	// initializes G_pDom, G_pMesh
	//loadSchemeEuler(g_genOpts.strScheme);
	//G_GhostMngEu.setDom(*G_pDomEu);

}