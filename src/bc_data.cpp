#include "bc_data.h"

#include "flow_params.h"

#include "IniFile.hpp"

t_BCList G_BCList;

// my_pvs - variables on the face of a real cell
// opp_pvs - variables on the face of a virtual cell 

// supersonic inflow bc (Dirichlet)
// IMPORTANT: prim vars here in global reference frame
void t_BCDataInflow::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {

	// rho
	opp_pvs[0] = 1.0;
	// u,v,w
	opp_pvs.setUVW(G_FreeStreamParams.getUInf());
	// p
	double M2 = std::pow(G_FreeStreamParams.getMach(), 2);
	double gMaMa = G_FlowModelParams.Gamma * M2;
	opp_pvs[4] = 1.0/gMaMa;

}

// supersonic outflow
// simple extrapolation, prim vars are in local or global rf
void t_BCDataOutFlow::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {

	opp_pvs = my_pvs;

}

// IMPORTANT: prim vars here in local reference frame (ksi, eta, dzeta)
// ksi is normal to the face
void t_BCDataEulerWall::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {

	opp_pvs = my_pvs;

	opp_pvs[1] *= -1.0;

}
// IMPORTANT: prim vars here in local reference frame (ksi, eta, dzeta)
// ksi is normal to the face
// Symmetry BC is identical to euler wall
void t_BCDataSym::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {

	opp_pvs = my_pvs;

	opp_pvs[1] *= -1.0;

}