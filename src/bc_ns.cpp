#include "bc_ns.h"

#include "flow_params.h"

#include "flux_ns.h"
#include "flow_model_perfect_gas.h"

#include <sstream>

//#pragma warning(push)
// Disable Visual Studio warnings
//#pragma warning(disable:4996)  // strtok may be unsafe ... ignoring)

t_BCListNS G_BCListNS;

// csv_my - variables on the face of a real cell
// csv_virt - variables on the face of a virtual cell 

// BC List

std::string t_BCListNS::getSupportedBCsStr() const {

	std::ostringstream ostr;

	ostr << t_BCNSInflowSup::getBCKindNameStat() << ", ";
	ostr << t_BCNSOutflowSup::getBCKindNameStat() << ", ";
	ostr << t_BCNSWall::getBCKindNameStat() << ", ";
	ostr << t_BCNSSym::getBCKindNameStat();

	return ostr.str();

}

t_FaceBCID t_BCListNS::getBCID(std::string fam_name) const {

	t_FaceBCID face_id;

	bool ok = false;

	for (int i = 0; i < _pBCs.size(); i++) {
		if (fam_name.compare(_pBCs[i]->getFamName()) == 0) {
			face_id.set(i);
			ok = true;
			break;
		}
	}
	if (!ok)
		hsLogError("t_BCListEuler:getBCID: failed to find bc with name=%s", fam_name.c_str());

	if (face_id.get() == t_FaceBCID::Fluid)
		hsLogError("BCListEuler: identificator of bc coincide with id of fluid face! Fix it!");

	return face_id;

};

const t_BCNSBase* t_BCListNS::getBC(int BCID) const {

#ifdef _DEBUG
	if (BCID<0 || BCID>_pBCs.size() - 1)
		hsLogError("t_BCListEuler::getBC: BCID=%d out of range", BCID);
#endif

	return _pBCs[BCID];

};

void t_BCListNS::addBCsetByName(std::string bc_set_name, std::string bc_kind_name, std::string& ini_data) {

	t_BCNSBase* pBC = nullptr;

	if (bc_kind_name.compare(t_BCNSInflowSup::getBCKindNameStat()) == 0)
		pBC = new t_BCNSInflowSup(bc_set_name);

	if (bc_kind_name.compare(t_BCNSOutflowSup::getBCKindNameStat()) == 0)
		pBC = new t_BCNSOutflowSup(bc_set_name);

	if (bc_kind_name.compare(t_BCNSWall::getBCKindNameStat()) == 0)
		pBC = new t_BCNSWall(bc_set_name);

	if (bc_kind_name.compare(t_BCNSSym::getBCKindNameStat()) == 0)
		pBC = new t_BCNSSym(bc_set_name);

	if (pBC == nullptr)
		hsLogMessage("Error:t_BCList:Unknown bc type!");
	else {
		hsLogMessage("Initializing bc set: %s", &(pBC->get_name()[0]));
		pBC->init(ini_data, "");
		_pBCs.push_back(pBC);
	}

}

bool t_BCListNS::getBCKindByFamName(const std::string& fam_name, t_BCKindNS& bc_kind) const {

	std::string bc_kind_str = "";
	bool ok = false;

	for (auto elem : _pBCs) {
		if (fam_name.compare(elem->getFamName()) == 0) {
			bc_kind_str = elem->getBCKindName();
		}
	};

	if (bc_kind_str.compare("") == 0) {
		// debug
		//hsLogMessage(
		//	"Warning:t_BCList::getBCKindBySectName: can't find bc for the section %s", &sect_name[0]);
	}
	else {
		ok = true;
		if (bc_kind_str.compare(t_BCNSInflowSup::getBCKindNameStat()) == 0) {
			bc_kind = t_BCKindNS::InflowSup;
		}
		if (bc_kind_str.compare(t_BCNSOutflowSup::getBCKindNameStat()) == 0) {
			bc_kind = t_BCKindNS::OutflowSup;
		}
		if (bc_kind_str.compare(t_BCNSWall::getBCKindNameStat()) == 0) {
			bc_kind = t_BCKindNS::WallNoSlip;
		}
		if (bc_kind_str.compare(t_BCNSSym::getBCKindNameStat()) == 0) {
			bc_kind = t_BCKindNS::Sym;
		}

	}
	return ok;

};

bool t_BCListNS::has(std::string sectionname) const {

	// unused
	t_BCKindNS bc_kind;

	return getBCKindByFamName(sectionname, bc_kind);

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

//#pragma warning(pop)