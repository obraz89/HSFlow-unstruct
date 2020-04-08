#pragma once

#include "flux_euler.h"

void calcWaveSpeedDavisEstim(
	const t_PrimVars& pvl, const t_PrimVars& pvr, double& sl, double& sr);

void calcRusanovFlux(
	const t_PrimVars& pvl, const t_PrimVars& pvr, t_Flux& flux);



