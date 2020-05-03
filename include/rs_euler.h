#pragma once

#include "flux_euler.h"

void calcWaveSpeedDavisEstim(
	const t_PrimVars& pvl, const t_PrimVars& pvr, double& sl, double& sr);

double calcWaveSpeedDavisEstimMaxAbs(const t_PrimVars& pvl, const t_PrimVars& pvr);

// calculate flux
void calcRSFlux(const t_PrimVars& pvl, const t_PrimVars& pvr, t_Flux& flux);


