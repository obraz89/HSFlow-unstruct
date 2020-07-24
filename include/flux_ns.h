#pragma once

#include "flow_vars_uvwpt.h"

#include "matrix_small.h"

//********************************Viscous Flux addition
// compute viscous flux vector
// all input vectors must be in global reference frame
// gradients matrix must be for UVWPT (NOT prim vec RUVWP)
void calcNSViscFlux(const t_Vec3& normal, const t_PrimVars& PV, const t_Mat<5, 3> GradUVWPT, t_VecConsVars& flux);

// calc CV & Total Flux for the face from primitive Vars 
// Total Flux = Euler Flux + Viscous Flux
// All input values in global reference frame
void calcNSFluxTotal(const t_PrimVars& pv, t_ConsVars& cv, t_VecConsVars& f);


