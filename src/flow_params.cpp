#include "flow_params.h"

t_FlowParamsFreeStream G_FreeStreamParams;

void initialize_flow_params() {

	G_FreeStreamParams.Mach = 2.0;

	G_FreeStreamParams.TinfDim = 300.0;

}