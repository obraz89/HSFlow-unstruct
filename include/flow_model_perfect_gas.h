#pragma once

#include "matrix_small.h"

struct t_FlowModelParams {

	double Gamma;

	double Pr;

};

extern t_FlowModelParams G_FlowModelParams;

void initialize_flow_model();



// non-dimensional speed of sound
double calcSoundSpeedByRP(double rho, double pressure); 

double calcTempByRP(double rho, double pressure);

double calcPressureByRT(double rho, double temperature);

double calcRhoByPT(double pressure, double temperature);

// gamma*Mach*Mach
double calcGMaMa();

// (gamma-1)*Mach*Mach
double calcG_1MaMa();

double calcViscosity(double T);

double calcThermalConductivity(double T);