#pragma once

#include "flow_vars_uvwpt.h"

#include "matrix_small.h"

//********************************Viscous Flux addition
// compute viscous flux vector
// all input vectors must be in global reference frame
// gradients matrix must be for UVWPT (NOT prim vec RUVWP)
void calcNSViscFlux(const t_Vec3& normal, const t_PrimVars& PV, const t_Mat<3, NConsVars> GradUVWPT, t_VecConsVars& flux);


