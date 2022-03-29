#include "schemes_euler.h"

#include "dom_euler_base.h"

#include "dom_euler_1st.h"

#include "dom_euler_lsq.h"

#include "dom_euler_1st_impl.h"

#include "dom_euler_lw.h"

static const char* strSupportedSchemes[] = {
	// first order explicit Godunov-type scheme
	"euler_1st",
	// second order explicit Godunov-type scheme with lsq reconstruction
	"euler_lsq",
	// first order implicit Godunov-type scheme
	"euler_1st_impl",
	// second order explicit Lax-Wendroff scheme
	"euler_lw" };

std::string getSupportedSchemesStr() {

	std::string s_all;

	for (auto ps : strSupportedSchemes) {
		s_all += std::string(ps)+", ";
	}

	return s_all;
}

void loadSchemeEuler(std::string strName) {

	// Euler 1st
	if (strName.compare(strSupportedSchemes[0]) == 0) {
		G_pDomEu = new t_DomEu1st();
		G_pDom = G_pDomEu;
		return;
	}

	// Euler lsq
	if (strName.compare(strSupportedSchemes[1]) == 0) {
		G_pDomEu = new t_DomEuLSQ();
		G_pDom = G_pDomEu;
		return;
	}

	// Euler 1st implicit
	if (strName.compare(strSupportedSchemes[2]) == 0) {
		G_pDomEu = new t_DomEu1stImpl();
		G_pDom = G_pDomEu;
		return;
	}

	// Euler Lax Wendroff
	if (strName.compare(strSupportedSchemes[3]) == 0) {
		G_pDomEu = new t_DomEuLW();
		G_pDom = G_pDomEu;
		return;
	}

	hsLogError("loadSchemeEuler: unsupported scheme=%s", strName.c_str());
	hsLogError("loadSchemeEuler: supported schemes are:%s", getSupportedSchemesStr().c_str());
	hsflow::TLog::flush();

};