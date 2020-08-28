#pragma once

#include "matrix_small.h"

#include "settings.h"

struct t_ViscType : public t_EnumStr {
	static const int Suther = 0;
	static const int Power = 1;
	void initValsStr() {
		ValsStr.push_back("Suther");
		ValsStr.push_back("Power");
	}

	bool operator==(int v) { return Val == v; }

	const std::string& defaultValStr() const { return ValsStr[0]; }

	t_ViscType() { initValsStr(); }
};

struct t_FlowModelParams {

	double Gamma;

	double Pr;

	t_ViscType ViscType;

};

extern t_FlowModelParams G_FlowModelParams;

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