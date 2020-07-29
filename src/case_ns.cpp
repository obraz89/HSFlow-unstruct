#include "settings.h"

#include "bc_ns.h"

#include "flow_params.h"

#include "flow_model_perfect_gas.h"

#include "schemes_ns.h"

#include "ghost_ns.h"

t_DomNSBase* G_pDomNS;

void load_case_ns(std::string& ini_data) {

	G_FreeStreamParams.init(ini_data, "");

	G_pBCList = &G_BCListNS;

	G_BCListNS.init(ini_data, "");

	G_pGhostMngBase = &G_GhostMngNS;
	// initializes G_pDom, G_pMesh
	loadSchemeNS(g_genOpts.strScheme);
	G_GhostMngNS.setDom(*G_pDomNS);

}