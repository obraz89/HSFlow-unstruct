#include "bc_euler.h"

#include "flow_params.h"

#include "flux_euler.h"
#include "flow_model_perfect_gas.h"

#include <sstream>

#pragma warning(push)
// Disable Visual Studio warnings
#pragma warning(disable:4996)  // strtok may be unsafe ... ignoring)

const std::string t_BCDataInflow::bc_kind = "bc_inflow";
const std::string t_BCDataOutFlow::bc_kind = "bc_outflow";
const std::string t_BCDataEulerWall::bc_kind = "bc_euler_wall";
const std::string t_BCDataSym::bc_kind = "bc_sym";

t_BCListEuler G_BCListEuler;

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

std::string t_BCListEuler::getSupportedBCsStr() {

	std::ostringstream ostr;

	ostr << t_BCDataInflow::bc_kind <<", ";
	ostr << t_BCDataOutFlow::bc_kind <<", ";
	ostr << t_BCDataEulerWall::bc_kind<<", ";
	ostr << t_BCDataSym::bc_kind;

	return ostr.str();

}

t_FaceBCID t_BCListEuler::getBCID(std::string sectionname) {
	t_BCKindEuler bc_kind;
	getBCKindBySectName(sectionname, bc_kind);

	t_FaceBCID faceBCID;

	if ((int)bc_kind == t_FaceBCID::Fluid)
		hsLogError("BCListEuler: identificator of bc coincide with id of fluid face! Fix it!");

	faceBCID.set((int)bc_kind);

	return faceBCID;

};

void t_BCListEuler::addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data) {

	t_BCDataFace* pBC = nullptr;

	if (bc_kind_name.compare(t_BCDataInflow::bc_kind) == 0)
		pBC = new t_BCDataInflow(bc_set_name);

	if (bc_kind_name.compare(t_BCDataOutFlow::bc_kind) == 0)
		pBC = new t_BCDataOutFlow(bc_set_name);

	if (bc_kind_name.compare(t_BCDataEulerWall::bc_kind) == 0)
		pBC = new t_BCDataEulerWall(bc_set_name);

	if (bc_kind_name.compare(t_BCDataSym::bc_kind) == 0)
		pBC = new t_BCDataSym(bc_set_name);

	if (pBC == nullptr)
		hsLogMessage("Error:t_BCList:Unknown bc type!");
	else {
		hsLogMessage("Initializing bc set: %s", &(pBC->get_name()[0]));
		pBC->init(ini_data, "");
		_pBCs.emplace(std::make_pair(bc_set_name, pBC));
	}

}

bool t_BCListEuler::getBCKindBySectName(const std::string& sect_name, t_BCKindEuler& bc_kind) {

	std::string bc_kind_str="";
	bool ok = false;

	for (auto elem : _pBCs) {
		if (sect_name.compare(elem.second->getSectName())==0) {
			bc_kind_str = elem.second->getBCKindName();
		}
	};

	if (bc_kind_str.compare("") == 0){
		// debug
		//hsLogMessage(
		//	"Warning:t_BCList::getBCKindBySectName: can't find bc for the section %s", &sect_name[0]);
	}
	else {
		ok = true;
		if (bc_kind_str.compare(t_BCDataInflow::bc_kind) == 0) {
			bc_kind = t_BCKindEuler::Inflow;
		}
		if (bc_kind_str.compare(t_BCDataOutFlow::bc_kind) == 0) {
			bc_kind = t_BCKindEuler::Outflow;
		}
		if (bc_kind_str.compare(t_BCDataEulerWall::bc_kind) == 0) {
			bc_kind = t_BCKindEuler::Wall;
		}
		if (bc_kind_str.compare(t_BCDataSym::bc_kind) == 0) {
			bc_kind = t_BCKindEuler::Sym;
		}

	}
	return ok;

};

bool t_BCListEuler::has(std::string sectionname) {

	// unused
	t_BCKindEuler bc_kind;
	
	return getBCKindBySectName(sectionname, bc_kind);

}

/*
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
*/

#pragma warning(pop)