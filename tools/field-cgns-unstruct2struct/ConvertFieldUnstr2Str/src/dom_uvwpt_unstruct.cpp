#include "mpi.h"

#include "dom_uvwpt_unstruct.h"

#include "io_field_unstruct.h"

#include "settings.h"

#include "cgnslib.h"
// TODO: make cgns wrapper as in HSFlow
#if CG_BUILD_SCOPE == 0
#undef  CG_BUILD_SCOPE
#define CG_BUILD_SCOPE 1
#endif

t_Dom5 G_DomUnst;

void t_Dom5::allocateFlowSolution() {

	ZonesSol = new t_ZoneFlowData[nZones];

	for (int i = iZneMPIs; i <= iZneMPIe; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowData& fdata = ZonesSol[i];

		lint nFaces = zne.getNFaces();
		lint nCellsReal = zne.getnCellsReal();

		fdata.PV = new t_PrimVars[nCellsReal];


	}

}

t_Dom5::~t_Dom5() {

		for (int i = iZneMPIs; i <= iZneMPIe; i++) {

			t_Zone& zne = Zones[i];
			t_ZoneFlowData& fdata = ZonesSol[i];

			if (fdata.PV!=nullptr) delete[] fdata.PV;

		}

		if (ZonesSol!=nullptr) delete[] ZonesSol;

}

void t_Dom5::initializeFlow() {

	hsLogMessage("* Reading initial field '%s'...",
		g_Settings.strFieldFnUstr.c_str());

	double time = loadField(g_Settings.strFieldFnUstr);


}

std::vector<std::string> t_Dom5::getFuncNamesIO() const {
	std::vector<std::string> fnames;
	fnames.push_back("VelocityX");
	fnames.push_back("VelocityY");
	fnames.push_back("VelocityZ");
	fnames.push_back("Pressure");
	fnames.push_back("Temperature");

	return fnames;
}

void t_Dom5::getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const {

	const t_Zone& Zne = Zones[zoneID];

	double* data = Vals.data();

	bool isCoord = false;

	// coords
	{
		int iXYZ = -1;

		for (int i = 0; i < 3; i++)
			if (name.compare(g_cgCoordNames[i]) == 0) {

				if (Vals.size() != Zne.getnVerts())
					hsLogError("t_DomainEuler::getDataAsArr: coord size mismatch");

				isCoord = true;
				iXYZ = i;
				break;
			}

		if (isCoord) {
			for (int iVert = 0; iVert < Zne.getnVerts(); iVert++)
				data[iVert] = Zne.getVert(iVert).xyz[iXYZ];
			//return;
		}
	}

	// funcs
	if (!isCoord) {
		if (Vals.size() != Zne.getnCellsReal())
			hsLogError("t_DomainEuler::getDataAsArr: flow field size mismatch");

		std::vector<std::string> func_names = getFuncNamesIO();

		int func_id = -1;

		for (int i = 0; i < func_names.size(); i++)
			if (name.compare(func_names[i]) == 0)
				func_id = i;

		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {
			const t_PrimVars& csv = getCellPV(zoneID, iCell);
			t_PrimVars pvio(csv);
			// velocityX
			if (func_id == 0) {
				data[iCell] = pvio.getU();
				continue;
			}
			// velocity Y
			if (func_id == 1) {
				data[iCell] = pvio.getV();
				continue;
			}
			// velocity Z
			if (func_id == 2) {
				data[iCell] = pvio.getW();
				continue;
			}
			// Pressure
			if (func_id == 3) {
				data[iCell] = pvio.getP();
				continue;
			}
			// Temperature
			if (func_id == 4) {
				data[iCell] = pvio.getT();
				continue;
			}

			hsLogError("t_DomainEuler::getDataAsArr: unknown name=%s", name.c_str());

		}

	}

	// debug
	//hsLogMessage("Array from dom_euler: name=%s", name.c_str());
	//for (int i = 0; i < Vals.size(); i++)
	//	std::cout << "Val #" << i << "=" << data[i] << ";";
	//std::cout << "\n";

}

double t_Dom5::loadField(std::string path_field) {

	double time = -1.0;

	//TLogSyncGuard logGuard;

	int f = -1;
	const int iBase = 1;  // assume only one base in the file

	char ok = 1;
	if (G_State.mpiRank == 0) do
	{
		ok = 0;

		if (cg_open(path_field.c_str(), CG_MODE_READ, &f) != CG_OK)
		{
			hsLogError("Can't open field file '%s' for reading: %s",
				path_field.c_str(), cg_get_error());
			break;
		}

		char szName[33];

		// Space dimensions
		int dimCell = 0, dimPhys = 0;
		if (cg_base_read(f, iBase, szName, &dimCell, &dimPhys) != CG_OK)
		{
			hsLogError("Can't read CGNS base node from '%s' ( %s )",
				path_field.c_str(), cg_get_error());
			break;
		}

		if (dimCell != g_pDomUnst->nDim) {
			hsLogError("CGNS: Inconsistent space dimensions");
			break;
		}

		// Number of zones (aka blocks)
		int nZones = 0;  cg_nzones(f, iBase, &nZones);
		if (nZones != g_pDomUnst->nZones) {
			hsLogError("CGNS: Inconsistent number of zones");
			break;
		}

		// Get time
		time = 0.0;
		do {
			//int nTmSteps = 0;  cg_biter_read(f, iBase, szName, &nTmSteps);
			//if (nTmSteps < 1) break;

			if (cg_goto(f, iBase, "BaseIterativeData_t", 1, NULL) != CG_OK)
				break;

			int nArrs = 0;  cg_narrays(&nArrs);
			if (nArrs < 1)  break;

			for (int ai = 1; ai <= nArrs; ++ai)
			{
				CG_DataType_t type;  int dim;  cgsize_t len[3];
				cg_array_info(ai, szName, &type, &dim, len);
				if (strcmp(szName, "TimeValues") == 0) {
					cg_array_read_as(ai, CG_RealDouble, &time);
					break;
				}
			}
		} while (false);

		G_State.time = time;

		ok = 1;
	} while (false); // if( G_State.mpiRank == 0 )

	MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, MPI_COMM_WORLD);
	if (!ok)
		return -1;

	MPI_Bcast(&time, 1, MPI_DOUBLE, 0/*root*/, MPI_COMM_WORLD);

	// Loop through zones
	for (int zi = 0; zi < nZones; ++zi)
	{
		t_Zone& zne = Zones[zi];
		const int NCellsReal = zne.getnCellsReal();

		// Data of the zone excluding ghosts
		TpakArraysDyn<double> newField;

		if ((this->iZneMPIs <= zi && zi <= this->iZneMPIe) || G_State.mpiRank == 0)
		{
			newField.reset(NConsVars, NCellsReal);
		}

		const int mpiTag = 'f' + 'l' + 'd' + zi;

		if (G_State.mpiRank == 0)
		{
			if (!read_zone_cgns(f, iBase, zi, newField))  ok = 0;

			// Send data from root to the zone's owner
			const int& rankDst = G_State.map_zone2rank[zi];
			if (rankDst != 0)  // don't send to myself
				MPI_Ssend(newField.data(), (ok ? newField.size() : 0), MPI_DOUBLE, rankDst, mpiTag, MPI_COMM_WORLD);
		}
		else if (g_pDomUnst->iZneMPIs <= zi && zi <= g_pDomUnst->iZneMPIe)
		{
			MPI_Recv(newField.data(), newField.size(), MPI_DOUBLE, 0/*root*/, mpiTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// TODO: Copy field to internal structs

		if (this->iZneMPIs <= zi && zi <= this->iZneMPIe) {

			for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

				t_PrimVars pvio;

				pvio.getU() = newField[0][iCell];
				pvio.getV() = newField[1][iCell];
				pvio.getW() = newField[2][iCell];
				pvio.getP() = newField[3][iCell];
				pvio.getT() = newField[4][iCell];

				getCellPV(zi, iCell) = pvio;

			}

		}



	}

	MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, MPI_COMM_WORLD);
	if (!ok)
		return -1;

	if (G_State.mpiRank == 0)
		cg_close(f);

	return time;

}

t_PrimVars t_Dom5::calcVertexPV(int iZone, int iVert) const{

		t_Zone& zne = Zones[iZone];

		t_PrimVars pv_vert;
		t_PrimVars pv_neig;

		pv_vert.reset();

		t_Vert& vert = zne.getVert(iVert);

		for (int j = 0; j < vert.NNeigCells; j++) {

				pv_neig = getCellPV(iZone, vert.pNeigCells[j]->Id);

				for (int k = 0; k < NConsVars; k++)
					pv_vert[k] += vert.pNeigCoefs[j] * pv_neig[k];

		}

		return pv_vert;

}