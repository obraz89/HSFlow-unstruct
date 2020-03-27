#include "common_data.h"

#include <cmath>

#include "flow_model.h"

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



	for (int i = 0; i < nCells; i++) {

		const t_Cell& cell = getCell(i);

		t_ConsVars dU;

		for (int j = 0; j < cell.NFaces; j++) {

			const t_Face& face = cell.getFace(j);

			double coef = cell.isMyFace(j) ? 1.0 : -1.0;

			coef *= dt*face.Area / cell.Volume;

			dU += coef*face.Flux;

		}

	}

}

void t_Zone::calcFaceFlux(lint iFace) {



}