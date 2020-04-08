#include "dom_euler.h"

#include "ghost_manager.h"

#include "rs_euler.h"

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

	// set real cell values
	for (int iZone = 0; iZone < nZones; iZone++) {

		t_Zone& zne = Zones[iZone];

		t_Cell* pcell;

		for (int i = 0; i < zne.getnCellsReal(); i++) {

			pcell = zne.getpCell(i);
			getCellCSV(iZone, i).setValAtInf();
		}

	}
	// set ghost values
	G_GhostManager.exchangeCSV();

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

	if (face.BCKind == t_FaceBCKind::Fluid) {

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

	if (face.BCKind == t_FaceBCKind::Inflow) {
		hsLogMessage("Inflow BC Flux: not implemented");
		return;
	}

	if (face.BCKind == t_FaceBCKind::Outflow) {
		hsLogMessage("OutFlow BC Flux: not implemented");
		return;
	}

	if (face.BCKind == t_FaceBCKind::Sym) {
		hsLogMessage("Sym BC Flux: not implemented");
		return;
	}

	if (face.BCKind == t_FaceBCKind::Wall) {
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