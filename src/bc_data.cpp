#include "bc_data.h"

#include "flow_params.h"

#include "IniFile.hpp"

#include <sstream>

#pragma warning(push)
// Disable Visual Studio warnings
#pragma warning(disable:4996)  // strtok may be unsafe ... ignoring)

const std::string t_BCDataInflow::bc_kind = "bc_inflow";
const std::string t_BCDataOutFlow::bc_kind = "bc_outflow";
const std::string t_BCDataEulerWall::bc_kind = "bc_euler_wall";
const std::string t_BCDataSym::bc_kind = "bc_sym";

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

// BC List

std::string t_BCList::getSupportedBCsStr() {

	std::ostringstream ostr;

	ostr << t_BCDataInflow::bc_kind <<", ";
	ostr << t_BCDataOutFlow::bc_kind <<", ";
	ostr << t_BCDataEulerWall::bc_kind<<", ";
	ostr << t_BCDataSym::bc_kind;

	return ostr.str();

}

std::vector<std::string> tokenize_str(std::string str) {

	std::vector<std::string> v;

	char* pch;
	pch = strtok(&str[0], " ,");
	while (pch != NULL)
	{
		v.push_back(pch);
		pch = strtok(NULL, " ,.-");
	}

	return v;
}

#pragma warning(pop)