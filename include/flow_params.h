#pragma once

struct t_FlowParamsFreeStream {

	double Mach;

	double TinfDim;

};

void initialize_flow_params();

extern t_FlowParamsFreeStream G_FreeStreamParams;
