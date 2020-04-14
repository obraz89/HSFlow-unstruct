#include "dom_euler.h"

#include "ghost_euler.h"

#include "rs_euler.h"

#include "bc_euler.h"

#include "settings.h"

#include "common_data.h"

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
	G_State.time =
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
	//hsLogMessage("Array from dom_euler: name=%s", name.c_str());
	//for (int i = 0; i < Vals.size(); i++)
	//	std::cout << "Val #" << i << "=" << data[i] << ";";
	//std::cout << "\n";

}

double t_DomainEuler::loadField(std::string path_field) {

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

	//MPI_Bcast(&ok, 1, MPI_CHAR, 0/*root*/, PETSC_COMM_WORLD);
	if (!ok)
		return -1;

	//MPI_Bcast(&time, 1, MPI_DOUBLE, 0/*root*/, PETSC_COMM_WORLD);

	// Loop through zones
	for (int zi = 0; zi < nZones; ++zi)
	{
		t_Zone& zne = Zones[zi];
		const int NCellsReal = zne.getnCellsReal();

		// Data of the zone excluding ghosts!!!
		TpakArraysDyn<double> newField;

		if ( G_State.mpiRank == 0)
		{
			newField.reset(NConsVars, NCellsReal);
		}
		if (G_State.mpiRank == 0)
		{
			if (!read_zone_cgns(f, iBase, zi, newField))  ok = 0;

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

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			t_PrimVarsIO pvio;
			
			pvio.u = newField[0][iCell];
			pvio.v = newField[1][iCell];
			pvio.w = newField[2][iCell];
			pvio.p = newField[3][iCell];
			pvio.t = newField[4][iCell];

			getCellCSV(zi, iCell) = pvio.calcConsVars();

		}

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

	double dU_max = 0.0;

	hsLogMessage("Computed dt=%lf", dt);

	G_State.ResidTot = 0.0;

	G_State.ResidNormVeloWall = 0.0;

	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		for (int iFace = 0; iFace < zne.getNFaces(); iFace++) {
			calcFaceFlux(iZone, iFace);
		}

		t_ConsVars dU;

		for (int iCell = 0; iCell < zne.getnCellsReal(); iCell++) {

			dU.reset();

			t_Cell& cell = zne.getCell(iCell);

			for (int j = 0; j < cell.NFaces; j++) {

				const t_Face& face = cell.getFace(j);

				// U[n+1] = U[n] - summ(flux_j), so we summ fluxes multiplied by minus 1
				double coef = cell.isMyFace(j) ? -1.0 : 1.0;

				coef *= dt * face.Area / cell.Volume;

				t_Flux flux = getFlux(iZone, face.Id);

				dU += coef * flux;

			}

			getCellCSV(iZone, iCell) += dU;

			double dU_norm = dU.norm();

			// max local resid for a cell
			if (dU_norm > dU_max) dU_max = dU_norm;
			// du/dt=rhs, sum all rhs to get resid
			G_State.ResidTot += dU_norm / dt;

		}

	}

	G_GhostMngEu.exchangeCSV();
	hsLogMessage("============");
	hsLogMessage("Time=%.6lf, Resid=%.6e, Local Max dU=%.6e", G_State.time, G_State.ResidTot, dU_max);
	hsLogMessage("Max normal velo resid at wall(sym):%.6e", G_State.ResidNormVeloWall);

	G_State.time += dt;

}

void t_DomainEuler::calcFaceFlux(int iZone, lint iFace) {

	t_Zone& zne = Zones[iZone];
	t_Face& face = zne.getFace(iFace);

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;

	// local vars for csvs, do not modify cell csv here
	t_ConsVars csv_my = getCellCSV(iZone, face.pMyCell->Id);

	t_Flux flux;

	if (face.BCId.get() == t_FaceBCID::Fluid) {

		t_ConsVars csv_op = getCellCSV(iZone, face.pOppCell->Id);

		t_PrimVars pvl = csv_my.calcPrimVars();
		t_PrimVars pvr = csv_op.calcPrimVars();

		// rotate everything to local rf
		R.set(mat_rot_coefs);

		// debug
		//hsLogMessage("Face #%d:", iFace);
		//hsLogMessage(R.to_str().c_str());

		pvl.rotate(R);
		pvr.rotate(R);

		calcRusanovFlux(pvl, pvr, flux);

		// rotate flux back
		R.set_inv(mat_rot_coefs);
		flux.rotate(R);

		// set flux for the face
		getFlux(iZone, iFace) = flux;

		return;


	}

	// TODO: now if not fluid, it's a BC

	t_BCKindEuler bc_kind = G_BCListEuler.getKind(face.BCId.get());

	if (bc_kind == t_BCKindEuler::Inflow) {

		t_ConsVars csv_face;
		csv_face.setValAtInf();

		R.set(mat_rot_coefs);
		csv_face.rotate(R);
		
		flux.calc(csv_face);

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;
		return;

	}

	if (bc_kind == t_BCKindEuler::Outflow) {

		t_ConsVars csv_face;

		csv_face = csv_my;

		R.set(mat_rot_coefs);

		csv_face.rotate(R);

		flux.calc(csv_face);

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;
		return;

	}

	if ((bc_kind == t_BCKindEuler::Wall) || (bc_kind == t_BCKindEuler::Sym)) {

		t_ConsVars csv_virt;

		R.set(mat_rot_coefs);

		csv_my.rotate(R);

		csv_virt = csv_my;
		// TODOL is this correct ? 
		csv_virt[1]*= -1.0;

		t_PrimVars pv_loc_my = csv_my.calcPrimVars();

		t_PrimVars pv_loc_virt = csv_virt.calcPrimVars();

		calcRusanovFlux(pv_loc_my, pv_loc_virt, flux);

		// DEBUG, get normal velocity (normal to the wall)
		// if local rf this is u_ksi
		// IMPORTANT: instead of u_ksi, we now check (rho*u_ksi)
		// it is easy, TODO: compute u_ksi via flux vars
		double v_norm = flux[0];
		if (v_norm > G_State.ResidNormVeloWall) G_State.ResidNormVeloWall = v_norm;

		// TODO: is this correct
		//flux[1] = 0.0;

		R.set_inv(mat_rot_coefs);

		flux.rotate(R);

		getFlux(iZone, iFace) = flux;

		return;

	}


	hsLogError(
		"t_DomainEuler::calcFaceFlux: unknow bc kind : Zone #%ld, face #%ld",
		iZone, iFace);



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

void t_DomainEuler::checkFlow() {

	double VolMin = HUGE_VAL;
	double VolMax = 0.0;
	double AreaMin = HUGE_VAL;
	double AreaMax = 0.0;
	double DiaMin = HUGE_VAL;
	double DiaMax = 0.0;

	for (int iZone = 0; iZone < nZones; iZone++) {

		const t_Zone& Zne = Zones[iZone];

		for (int iCell = 0; iCell < Zne.getnCellsReal(); iCell++) {

			const t_Cell& cell = Zne.getCell(iCell);

			if (cell.Volume < VolMin) VolMin = cell.Volume;
			if (cell.Volume > VolMax) VolMax = cell.Volume;

			if (cell.Diameter < DiaMin) DiaMin = cell.Diameter;
			if (cell.Diameter > DiaMax) DiaMax = cell.Diameter;
		}

		for (int iFace = 0; iFace < Zne.getNFaces(); iFace++) {

			const t_Face& face = Zne.getFace(iFace);

			if (face.Area < AreaMin) AreaMin = face.Area;
			if (face.Area > AreaMax) AreaMax = face.Area;
		}
	}

	hsLogMessage("Domain info:");
	hsLogMessage("Min Volume=%lf, Max Volume=%lf", VolMin, VolMax);
	hsLogMessage("Min Cell Diameter=%lf, Max Cell Diameter=%lf", DiaMin, DiaMax);
	hsLogMessage("Min Face Area=%lf, Max Face Area=%lf", AreaMin, AreaMax);

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