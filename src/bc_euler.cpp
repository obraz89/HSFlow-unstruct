#include "bc_euler.h"

#include "flow_params.h"

#include "flux_euler.h"
#include "flow_model_perfect_gas.h"

#include <sstream>

#pragma warning(push)
// Disable Visual Studio warnings
#pragma warning(disable:4996)  // strtok may be unsafe ... ignoring)

const std::string t_BCDataInflow::bc_kind_str = "bc_inflow";
const t_BCKindEuler t_BCDataInflow::bc_kind = t_BCKindEuler::Inflow;

const std::string t_BCDataOutFlow::bc_kind_str = "bc_outflow";
const t_BCKindEuler t_BCDataOutFlow::bc_kind = t_BCKindEuler::Outflow;

const std::string t_BCDataEulerWall::bc_kind_str = "bc_euler_wall";
const t_BCKindEuler t_BCDataEulerWall::bc_kind = t_BCKindEuler::Wall;

const std::string t_BCDataSym::bc_kind_str = "bc_sym";
const t_BCKindEuler t_BCDataSym::bc_kind = t_BCKindEuler::Sym;

t_BCListEuler G_BCListEuler;

// csv_my - variables on the face of a real cell
// csv_virt - variables on the face of a virtual cell 

// supersonic inflow bc (Dirichlet)
// IMPORTANT: vars here in global reference frame
void t_BCDataInflow::yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const{

	t_PrimVars pv_virt;

	// rho
	pv_virt.setR(1.0);
	// u,v,w
	pv_virt.setUVW(G_FreeStreamParams.getUInf());
	// p
	pv_virt.setP(calcPressureByRT(1.0,1.0)) ;

	csv_virt = pv_virt.calcConsVars();

}

// supersonic outflow
// simple extrapolation, prim vars are in local or global rf
void t_BCDataOutFlow::yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const{

	csv_virt = csv_my;

}

// IMPORTANT: prim vars here in local reference frame (ksi, eta, dzeta)
// ksi is normal to the face
void t_BCDataEulerWall::yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const{

	csv_virt = csv_my;

	// TODO: this should tend to zero as field converges, check it
	csv_virt[1] *= -1.0;

}
// IMPORTANT: prim vars here in local reference frame (ksi, eta, dzeta)
// ksi is normal to the face
// Symmetry BC is identical to euler wall
void t_BCDataSym::yield(const t_ConsVars& csv_my, t_ConsVars& csv_virt) const{

	csv_virt = csv_my;

	csv_virt[1] *= -1.0;

}

// BC List

std::string t_BCListEuler::getSupportedBCsStr() const{

	std::ostringstream ostr;

	ostr << t_BCDataInflow::bc_kind_str <<", ";
	ostr << t_BCDataOutFlow::bc_kind_str <<", ";
	ostr << t_BCDataEulerWall::bc_kind_str <<", ";
	ostr << t_BCDataSym::bc_kind_str;

	return ostr.str();

}

t_FaceBCID t_BCListEuler::getBCID(std::string sectionname) const{

	t_FaceBCID face_id;

	bool ok = false;

	for (int i = 0; i < _pBCs.size(); i++) {
		if (sectionname.compare(_pBCs[i]->getSectName()) == 0) {
			face_id.set(i);
			ok = true;
			break;
		}
	}
	if (!ok)
		hsLogError("t_BCListEuler:getBCID: failed to find bc with name=%s", sectionname.c_str());

	if (face_id.get() == t_FaceBCID::Fluid)
		hsLogError("BCListEuler: identificator of bc coincide with id of fluid face! Fix it!");

	return face_id;

};

const t_BCEulerBase* t_BCListEuler::getBC(int BCID) const{

#ifdef _DEBUG
	if (BCID<0 || BCID>_pBCs.size() - 1)
		hsLogError("t_BCListEuler::getBC: BCID=%d out of range", BCID);
#endif

	return _pBCs[BCID];

};

void t_BCListEuler::addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data) {

	t_BCEulerBase* pBC = nullptr;

	if (bc_kind_name.compare(t_BCDataInflow::bc_kind_str) == 0)
		pBC = new t_BCDataInflow(bc_set_name);

	if (bc_kind_name.compare(t_BCDataOutFlow::bc_kind_str) == 0)
		pBC = new t_BCDataOutFlow(bc_set_name);

	if (bc_kind_name.compare(t_BCDataEulerWall::bc_kind_str) == 0)
		pBC = new t_BCDataEulerWall(bc_set_name);

	if (bc_kind_name.compare(t_BCDataSym::bc_kind_str) == 0)
		pBC = new t_BCDataSym(bc_set_name);

	if (pBC == nullptr)
		hsLogMessage("Error:t_BCList:Unknown bc type!");
	else {
		hsLogMessage("Initializing bc set: %s", &(pBC->get_name()[0]));
		pBC->init(ini_data, "");
		_pBCs.push_back(pBC);
	}

}

bool t_BCListEuler::getBCKindBySectName(const std::string& sect_name, t_BCKindEuler& bc_kind) const{

	std::string bc_kind_str="";
	bool ok = false;

	for (auto elem : _pBCs) {
		if (sect_name.compare(elem->getSectName())==0) {
			bc_kind_str = elem->getBCKindName();
		}
	};

	if (bc_kind_str.compare("") == 0){
		// debug
		//hsLogMessage(
		//	"Warning:t_BCList::getBCKindBySectName: can't find bc for the section %s", &sect_name[0]);
	}
	else {
		ok = true;
		if (bc_kind_str.compare(t_BCDataInflow::bc_kind_str) == 0) {
			bc_kind = t_BCKindEuler::Inflow;
		}
		if (bc_kind_str.compare(t_BCDataOutFlow::bc_kind_str) == 0) {
			bc_kind = t_BCKindEuler::Outflow;
		}
		if (bc_kind_str.compare(t_BCDataEulerWall::bc_kind_str) == 0) {
			bc_kind = t_BCKindEuler::Wall;
		}
		if (bc_kind_str.compare(t_BCDataSym::bc_kind_str) == 0) {
			bc_kind = t_BCKindEuler::Sym;
		}

	}
	return ok;

};

bool t_BCListEuler::has(std::string sectionname) const{

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