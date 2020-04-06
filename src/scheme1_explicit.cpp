#include "common_data.h"

#include <cmath>

#include "flow_model.h"

#include "rs_procs.h"

void t_Domain::makeTimeStep() {

	double dt=HUGE_VAL;

	for (int iZone = 0; iZone < nZones; iZone++) {

		double dt_z = Zones[iZone].calcDt();

		if (dt_z < dt) dt = dt_z;
	}

	for (int iZone = 0; iZone < nZones; iZone++)
		Zones[iZone].makeTimeStep(dt);

}

void t_Zone::makeTimeStep(double dt) {

	for (int i = 0; i < nFaces; i++) {
		calcFaceFlux(i);
	}

	t_ConsVars dU;

	for (int i = 0; i < nCellsReal; i++) {

		const t_Cell& cell = getCell(i);

		for (int j = 0; j < cell.NFaces; j++) {

			const t_Face& face = cell.getFace(j);

			double coef = cell.isMyFace(j) ? 1.0 : -1.0;

			coef *= dt*face.Area / cell.Volume;

			dU += coef*face.Flux;

		}

	}

}

void t_Zone::calcFaceFlux(lint iFace) {


	t_Face& face = Faces[iFace];

	t_MatRotN mat_rot_coefs;
	mat_rot_coefs.calc_rot_angles_by_N(face.Normal);

	t_SqMat3 R;
	R.set(mat_rot_coefs);

	hsLogMessage("Face #%d:", iFace);
	hsLogMessage(&R.to_str()[0]);

	if (face.BCKind == t_FaceBCKind::Fluid) {

		const t_PrimVars& pvl = face.pMyCell->PrimVars;
		const t_PrimVars& pvr = face.pOppCell->PrimVars;

		// rotate everything to local rf

		t_PrimVars pvl_loc = pvl;
		pvl_loc.rotate(R);

		t_PrimVars pvr_loc = pvr;
		pvr_loc.rotate(R);

		t_Flux flux_loc;

		calcRusanovFlux(pvl_loc, pvr_loc, flux_loc);

		// rotate flux back

		R.set_inv(mat_rot_coefs);

		face.Flux = flux_loc.rotate(R);

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

double t_Zone::calcDt() {

	return 0.0;

}