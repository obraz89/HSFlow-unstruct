#include "flow_common.h"

#include "flow_model_perfect_gas.h"

#include "flow_params.h"

t_FlowModelParams G_FlowModelParams;

void initialize_flow_model() {

	G_FlowModelParams.Gamma = 1.4;

	G_FlowModelParams.Pr = 0.72;

}



// common functions

double calcSoundSpeedByRP(double rho, double pressure) {

	return sqrt(G_FlowModelParams.Gamma * pressure / rho);

}

double calcGMaMaInv() {
	double M2 = G_FreeStreamParams.getMach() * G_FreeStreamParams.getMach();
	return 1.0 / (G_FlowModelParams.Gamma * M2);
}
