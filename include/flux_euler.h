#pragma once

#include "flow_vars_uvwpt.h"

class t_FluxEu : public t_VecConsVars{
public:
	t_FluxEu() :t_VecConsVars() {}
	t_FluxEu(const t_VecConsVars& v) :t_VecConsVars(v) {}
	void calc(const t_ConsVars& csv);
	void calc(const t_PrimVars& pvs);

};

// calc CV & Flux from primitive Vars
void calcCVFlux(const t_PrimVars& pv, t_ConsVars& cv, t_FluxEu& f);


