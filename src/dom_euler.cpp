#include "dom_euler.h"

#include "ghost_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "settings.h"

#include <fstream>

#include "io-field.h"


t_DomainEuler G_Domain;

void t_DomainEuler::allocateFlowSolution() {

	ZonesSol = new t_ZoneFlowData[nZones];

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowData& fdata = ZonesSol[i];
		
		lint nFaces = zne.getNFaces();
		lint nCellsTot = zne.getnCellsTot();

		fdata.Fluxes = new t_Flux[nFaces];
		fdata.ConsVars = new t_ConsVars[nCellsTot];


	}

}

void t_DomainEuler::initializeFlow() {

	double time = -1;

	if (g_genOpts.strInitFieldFN.empty()) {
		// set real cell values
		for (int iZone = 0; iZone < nZones; iZone++) {

			t_Zone& zne = Zones[iZone];

			t_Cell* pcell;

			for (int i = 0; i < zne.getnCellsReal(); i++) {

				pcell = zne.getpCell(i);
				getCellCSV(iZone, i).setValAtInf();
			}

		}

		time = 0.0;
	}
	else {

		hsLogMessage("* Reading initial field '%s'...",
			g_genOpts.strInitFieldFN.c_str());

		time = loadField(g_genOpts.strInitFieldFN);
		if (time < 0)
			hsLogError("Failed to load cgns field");

	}

	// Set initial time in time-stepping
	G_State.timeStart =
		(g_genOpts.timeStart >= 0) ? g_genOpts.timeStart : time;

	
	// set ghost values
	G_GhostMngEu.exchangeCSV();

}

std::vector<std::string> t_DomainEuler::getFuncNamesIO() {
	std::vector<std::string> fnames;
	fnames.push_back("VelocityX");
	fnames.push_back("VelocityY");
	fnames.push_back("VelocityZ");
	fnames.push_back("Pressure");
	fnames.push_back("Temperature");

	return fnames;
}

double t_DomainEuler::loadField(std::string fileName) {

	double time = -1.0;

	//TLogSyncGuard logGuard;

	int f = -1;
	const int iBase = 1;  // assume only one base in the file

	char ok = 1;
	if (G_State.mpiRank == 0) do
	{
		ok = 0;

		if (cg_open(fileName.c_str(), CG_MODE_READ, &f) != CG_OK)
		{
			hsLogError("Can't open field file '%s' for reading: %s",
				fileName.c_str(), cg_get_error());
			break;
		}

		char szName[33];

		// Space dimensions
		int dimCell = 0, dimPhys = 0;
		if (cg_base_read(f, iBase, szName, &dimCell, &dimPhys) != CG_OK)
		{
			hsLogError("Can't read CGNS base node from '%s' ( %s )",
				fileName.c_str(), cg_get_error());
			break;
		}

		if (dimCell != G_Domain.nDim) {
			hsLogError("CGNS: Inconsistent space dimensions");
			break;
		}

		// Number of zones (aka blocks)
		int nZones = 0;  cg_nzones(f, iBase, &nZones);
		if (nZones != G_Domain.nZones) {
			hsLogError("CGNS: Inconsistent number of zones");
			break;
		}

		// Get time
		time = 0.0;
		do {
			int nTmSteps = 0;  cg_biter_read(f, iBase, szName, &nTmSteps);
			if (nTmSteps < 1) break;

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

		ok = 1;
	} while (false); // if( G_State.mpiRank == 0 )

	//MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, PETSC_COMM_WORLD);
	if (!ok)
		return -1;

	//MPI_Bcast(&time, 1, MPI_DOUBLE, 0/*root*/, PETSC_COMM_WORLD);

	// Loop through zones
	for (int zi = 0; zi < nZones; ++zi)
	{
		t_Zone& zne = Zones[zi];
		const int NCellsReal = zne.getnCellsReal();
		const int NVerts = zne.getnVerts();

		// Data of the zone excluding ghosts!!!
		TpakArraysDyn<double> newField, newGrid;

		if ( G_State.mpiRank == 0)
		{
			newField.reset(NConsVars, NCellsReal);
			newGrid.reset(this->nDim, NVerts);
		}
		if (G_State.mpiRank == 0)
		{
			if (!read_zone_cgns(f, iBase, zi, newGrid, newField))  ok = 0;

			// Send data from root to the zone's owner
			//const int& rankDst = G_State.map_zone2rank[zi];
			//if (rankDst != 0)  // don't send to myself
			//	MPI_Ssend(newField.data(), (ok ? newField.size() : 0), MPI_DOUBLE, rankDst, mpiTag, PETSC_COMM_WORLD);
		}
		else if (G_Domain.bs <= zi && zi <= G_Domain.be)
		{
			//MPI_Recv(newField.data(), newField.size(), MPI_DOUBLE, 0/*root*/, mpiTag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// TODO: Copy field to internal structs

	}

	//MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, PETSC_COMM_WORLD);
	if (!ok)
		return -1;

	if (G_State.mpiRank == 0)
		cg_close(f);

	return time;

}

void t_DomainEuler::makeTimeStep() {

	double dt = calcDt();

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		for (int iFace = 0; iFace < zne.getNFaces(); iFace++) {
			calcFaceFlux(iZone, iFace);
		}

		t_ConsVars dU;

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			t_Cell& cell = zne.getCell(iCell);

			for (int j = 0; j < cell.NFaces; j++) {

				const t_Face& face = cell.getFace(j);

				double coef = cell.isMyFace(j) ? 1.0 : -1.0;

				coef *= dt * face.Area / cell.Volume;

				t_Flux flux = getFlux(iZone, face.Id);

				dU += coef * flux;

			}

			getCellCSV(iZone, iCell) += dU;

		}

	}

}

void t_DomainEuler::calcFaceFlux(int iZone, lint iFace) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;
	R.set(mat_rot_coefs);

	hsLogMessage("Face #%d:", iFace);
	hsLogMessage(&R.to_str()[0]);

	if (face.BCId.get() == t_FaceBCID::Fluid) {

		//t_PrimVars pvl = face.pMyCell->ConsVars.calcPrimVars();
		t_PrimVars pvl = getCellCSV(iZone, face.pMyCell->Id).calcPrimVars();
		//t_PrimVars pvr = face.pOppCell->ConsVars.calcPrimVars();
		t_PrimVars pvr = getCellCSV(iZone, face.pOppCell->Id).calcPrimVars();

		// rotate everything to local rf
		pvl.rotate(R);
		pvr.rotate(R);

		t_Flux flux_loc;

		calcRusanovFlux(pvl, pvr, flux_loc);

		// rotate flux back
		R.set_inv(mat_rot_coefs);
		flux_loc.rotate(R);
		//face.Flux = flux_loc.rotate(R);
		getFlux(iZone, iFace) = flux_loc;

		return;
	}

	if (face.BCId == (int)t_BCKindEuler::Inflow) {
		hsLogMessage("Inflow BC Flux: not implemented");
		return;
	}

	if (face.BCId == (int)t_BCKindEuler::Outflow) {
		hsLogMessage("OutFlow BC Flux: not implemented");
		return;
	}

	if (face.BCId == (int)t_BCKindEuler::Sym) {
		hsLogMessage("Sym BC Flux: not implemented");
		return;
	}

	if (face.BCId == (int)t_BCKindEuler::Wall) {
		hsLogMessage("Wall BC Flux: not implemented");
		return;
	}

	hsLogMessage("calcFaceFlux: unknown face type");
	return;


}

double t_DomainEuler::calcDt() {

	return 0.0;

}

void t_DomainEuler::dump_flow() {

	std::string fn("dump_flow.txt");

	std::ofstream ofstr(fn);

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		ofstr << "=========Zone #" << iZone << "===========\n";

		for (int i = 0; i < zne.getnCellsTot(); i++) {

			pcell = zne.getpCell(i);
			ofstr << "cell #" << i << getCellCSV(iZone, i).to_str();
		}

	}

	ofstr.flush();

}

void t_DomainEuler::dump_geom() {



}

t_DomainEuler::~t_DomainEuler() {

	for (int i = 0; i < nZones; i++) {

		t_Zone& zne = Zones[i];
		t_ZoneFlowData& fdata = ZonesSol[i];

		lint nFaces = zne.getNFaces();
		lint nCellsTot = zne.getnCellsTot();

		delete[] fdata.Fluxes;
		delete[] fdata.ConsVars;

	}

	delete[] ZonesSol;

}