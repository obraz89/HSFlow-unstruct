#pragma once

#include "dom_ns_base.h"

#include "dom_ns_lsq_expl.h"

static const char* strSupportedSchemes[] = { "ns_lsq"};

static void loadSchemeNS(std::string strName) {

	if (strName.compare(strSupportedSchemes[0]) == 0) {
		G_pDomNS = new t_DomNSLSQ();
		G_pDom = G_pDomNS;
		return;
	}

	hsLogError("loadSchemeNS: unsupported scheme=%s", strName.c_str());
	hsflow::TLog::flush();

};


