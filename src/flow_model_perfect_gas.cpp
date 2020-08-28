#include "flow_model_perfect_gas.h"

#include "flow_params.h"

t_FlowModelParams G_FlowModelParams;

// common functions

double calcSoundSpeedByRP(double rho, double pressure) {

	return sqrt(G_FlowModelParams.Gamma * pressure / rho);

}

double calcTempByRP(double rho, double pressure) {

	return calcGMaMa() * pressure / rho;

}

double calcPressureByRT(double rho, double T) {
	return 1.0 / calcGMaMa() * rho * T;
}

double calcRhoByPT(double pressure, double temperature) {
	return calcGMaMa() * pressure / temperature;
};

double calcGMaMa() {
	double M2 = G_FreeStreamParams.getMach() * G_FreeStreamParams.getMach();
	return G_FlowModelParams.Gamma * M2;
}

double calcG_1MaMa() {
	double M2 = G_FreeStreamParams.getMach() * G_FreeStreamParams.getMach();
	return (G_FlowModelParams.Gamma - 1) * M2;
}

// non-dimensional viscosity
double calcViscosityPower(double T) {
	return pow(T, 0.75);
};

// non-dimensional viscosity
double calcViscositySuther(double T) {

	const double _S = 110.4 / G_FreeStreamParams.getTinfDim();

	return (1 + _S) / (T + _S) * sqrt(T*T*T);

};

double calcViscosity(double T) {
	if (G_FlowModelParams.ViscType==t_ViscType::Power)
		return calcViscosityPower(T);
	if (G_FlowModelParams.ViscType==t_ViscType::Suther)
		return calcViscositySuther(T);
	hsLogError("CalcViscosity::unsupported visc calc type");
	return 0;
}

double calcThermalConductivity(double T) {

	return calcViscosity(T) / (calcG_1MaMa()*G_FlowModelParams.Pr);

};


