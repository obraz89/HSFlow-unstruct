#pragma once

#include "flux_euler.h"

class t_FluxNS : public t_VecConsVars {
public:
	t_FluxNS() :t_VecConsVars() {}
	t_FluxNS(const t_VecConsVars& v) :t_VecConsVars(v) {}
	void calc(const t_ConsVars& csv);
	void calc(const t_PrimVars& pvs);

};



// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_FluxNS& f);
