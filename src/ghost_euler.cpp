#include "mpi.h"

#include "ghost_euler.h"

t_GhostMngEuler G_GhostMngEu;

void  t_GhostMngEuler::exchangeCSV() {

	const int rankMy = G_State.mpiRank;

	for (int i = 0; i < _pDom->nZones; i++) {

		const int rankI = G_State.map_zone2rank[i];

		for (int j = 0; j < _pDom->nZones; j++) {

			const int rankJ = G_State.map_zone2rank[j];

			// don't send to myself
			if (i == j) continue;
			// not my pair, skip
			if (!(rankI == rankMy || rankJ == rankMy)) continue;

			// ghost layer for zone i from zone j
			t_GhostLayer* glayer = _pGLayers[getPlainInd(i, j)];

			// local exchange, rankI==rankJ==rankMy
			if (rankI == rankJ) {

				for (int k = 0; k < glayer->size(); k++) {

					lint offset = calcIndOffset(i, j);
					_pDomEu->getCellCSV(i, offset + k) = _pDomEu->getCellCSV(j, glayer->data[k].id_dnr);

				}

			}
			// mpi exchange
			else {

				const int mpiTag = 'c' + 's' + 'v' + getPlainInd(i, j);

				t_ArrDbl arr;
				int arr_size = NConsVars * glayer->size();
				arr.alloc(arr_size);

				if (rankMy == rankJ) {
					// prepare & send data for zone i from zone j

					for (int k = 0; k < glayer->size(); k++) {

						const t_ConsVars& csv = _pDomEu->getCellCSV(j, glayer->data[k].id_dnr);

						lint shift = NConsVars * lint(k);

						for (int l = 0; l < NConsVars; l++)
							arr.set(shift + l, csv[l]);

					}

					MPI_Ssend(arr.data(), arr_size, MPI_DOUBLE, rankI, mpiTag,
						MPI_COMM_WORLD);

				}
				else {	// rankMy == rankI
					// receive data from zone j

					MPI_Recv(arr.data(), arr_size, MPI_DOUBLE,
						rankJ, mpiTag, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

					lint offset = calcIndOffset(i, j);
					for (int k = 0; k < glayer->size(); k++) {

						t_ConsVars& csv = _pDomEu->getCellCSV(i, offset + k);

						lint shift = NConsVars * lint(k);

						for (int l = 0; l < NConsVars; l++)
							csv[l] = arr.data()[shift + l];

					}
				}

			}

		}
	}

	_pDomEu->checkMinMaxCSV();

};

void t_GhostMngEuler::exchangeReconstData() {
	hsLogMsgAllRanks("t_GhostMngEuler::exchangeReconstData not implemented yet, ghosts not receiving reconstr data");
}