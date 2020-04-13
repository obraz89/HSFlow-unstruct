#include "dom_euler.h"

#include "ghost_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "settings.h"

#include <fstream>


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

std::vector<std::string> t_DomainEuler::getFuncNamesIO() const{
	std::vector<std::string> fnames;
	fnames.push_back("VelocityX");
	fnames.push_back("VelocityY");
	fnames.push_back("VelocityZ");
	fnames.push_back("Pressure");
	fnames.push_back("Temperature");

	return fnames;
}

void t_DomainEuler::getDataAsArr(std::string name, int zoneID, t_ArrDbl& Vals) const {

	const t_Zone& Zne = Zones[zoneID];

	double* data = Vals.data();

	bool isCoord = false;

	// coords
	{
		 int iXYZ = -1;

		for (int i = 0; i < 3; i++)
			if (name.compare(g_cgCoordNames[i])==0) {

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
	if (!isCoord){
		if (Vals.size() != Zne.getnCellsReal())
			hsLogError("t_DomainEuler::getDataAsArr: flow field size mismatch");

		std::vector<std::string> func_names = getFuncNamesIO();

		int func_id = -1;

		for (int i = 0; i < func_names.size(); i++)
			if (name.compare(func_names[i]) == 0)
				func_id = i;

		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {
			const t_ConsVars& csv = getCellCSV(zoneID, iCell);
			t_PrimVarsIO pvio(csv);
			// velocityX
			if (func_id==0) {
				data[iCell] = pvio.u;
				continue;
			}
			// velocity Y
			if (func_id==1) {
				data[iCell] = pvio.v;
				continue;
			}
			// velocity Z
			if (func_id==2) {
				data[iCell] = pvio.w;
				continue;
			}
			// Pressure
			if (func_id==3) {
				data[iCell] = pvio.p;
				continue;
			}
			// Pressure
			if (func_id==4) {
				data[iCell] = pvio.t;
				continue;
			}

			hsLogError("t_DomainEuler::getDataAsArr: unknown name=%s", name.c_str());

		}
		
	}

	// debug
	hsLogMessage("Array from dom_euler: name=%s", name.c_str());
	for (int i = 0; i < Vals.size(); i++)
		std::cout << "Val #" << i << "=" << data[i] << ";";
	std::cout << "\n";

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

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);
	t_ConsVars csv_op; 

	// TODO: assuming that if not fluid face, then it is a bc face

	if (face.BCId.get() == t_FaceBCID::Fluid) {

		csv_op = getCellCSV(iZone, face.pOppCell->Id);
	}
	else do {

		t_BCKindEuler bc_kind = G_BCListEuler.getKind(face.BCId.get());


		if ((bc_kind == t_BCKindEuler::Inflow) ||
			(bc_kind == t_BCKindEuler::Outflow)) {
			// not rotating, setting bc in glob rf
			G_BCListEuler.getBC(face.BCId.get())->yield(csv_my, csv_op);
			break;
		}

		if ((bc_kind == t_BCKindEuler::Sym) ||
			(bc_kind == t_BCKindEuler::Wall)) {
			// rotate
			R.set(mat_rot_coefs);
			csv_my.rotate(R);
			csv_op.rotate(R);
			// set BC in local reference frame
			G_BCListEuler.getBC(face.BCId.get())->yield(csv_my, csv_op);
			// rotate back
			R.set_inv(mat_rot_coefs);
			csv_my.rotate(R);
			csv_op.rotate(R);

			break;

		}

		hsLogError(
			"t_DomainEuler::calcFaceFlux: unknow bc kind : Zone #%ld, face #%ld",
			iZone, iFace);
	} while (false);

	// My CSV & Opp CSV are set, calculate flux for fluid face
	t_PrimVars pvl = csv_my.calcPrimVars();
	t_PrimVars pvr = csv_op.calcPrimVars();

	// rotate everything to local rf
	R.set(mat_rot_coefs);

	hsLogMessage("Face #%d:", iFace);
	hsLogMessage(R.to_str().c_str());

	pvl.rotate(R);
	pvr.rotate(R);

	t_Flux flux_loc;

	calcRusanovFlux(pvl, pvr, flux_loc);

	// rotate flux back
	R.set_inv(mat_rot_coefs);
	flux_loc.rotate(R);
	
	// set flux for the face
	getFlux(iZone, iFace) = flux_loc;

	return;


}

double t_DomainEuler::calcDt() const{

	double dt = HUGE_VAL;

	for (int iZone = 0; iZone < nZones; iZone++) {

		const t_Zone& Zne = Zones[iZone];

		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {

			const t_Cell& cell = Zne.getCell(iCell);

			const t_ConsVars& csv = getCellCSV(iZone, iCell);

			t_PrimVars pv = csv.calcPrimVars();

			double c = calcSoundSpeedByRP(pv.getR(), pv.getP());

			// rough estimates of |U-c|, |U| and |U+c| is just abs(U)+c
			double v_max = pv.getUVW().norm() + c;

			double dt_cur = g_genOpts.CFL * cell.Diameter / v_max;

			if (dt_cur < dt) dt = dt_cur;

		}

	}

	return dt;

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