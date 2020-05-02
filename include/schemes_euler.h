#pragma once

#include "dom_euler_base.h"

#include "dom_euler_1st.h"

#include "dom_euler_lsq.h"

#include "dom_euler_1st_impl.h"

static const char* strSupportedSchemes[] = { "euler_1st", "euler_lsq", "euler_1st_impl" };

static void loadSchemeEuler(std::string strName) {

	if (strName.compare(strSupportedSchemes[0]) == 0) {
		G_pDom = new t_DomEu1st();
		G_pMesh = G_pDom;
		return;
	}

	if (strName.compare(strSupportedSchemes[1]) == 0) {
		G_pDom = new t_DomEuLSQ();
		G_pMesh = G_pDom;
		return;
	}

	if (strName.compare(strSupportedSchemes[2]) == 0) {
		G_pDom = new t_DomEu1stImpl();
		G_pMesh = G_pDom;
		return;
	}

	hsLogError("loadSchemeEuler: unsupported scheme=%s", strName.c_str());
	hsflow::TLog::flush();

};

