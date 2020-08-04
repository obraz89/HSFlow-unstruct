#include "common_data_struct.h"

#include "logging.h"

#include "mpi.h"

bool assignZonesToProcs()
{
	if (G_Domain.nZones < G_State.mpiNProcs)
	{
		hsLogError(
			"Domain zones count (%d) is less than MPI procs (%d)",
			G_Domain.nZones, G_State.mpiNProcs
		);
		return false;
	}

	G_State.map_zone2rank = new int[G_Domain.nZones];
	for (int b = 0; b < G_Domain.nZones; ++b)  G_State.map_zone2rank[b] = -1;

	G_Domain.map_iZne2cgID = new int[G_Domain.nZones];

	hsLogMessage(" ");
	hsLogMessage("* Grid zones distribution through procs");

	if (G_State.mpiRank == 0)
		do
		{
			hsLogMessage(" struct: default zones per cpu layout");

			//
			// Default layout: one-to-one mapping of zones indices to CGNS zone IDs
			//
			for (int b = 0; b < G_Domain.nZones; ++b)
				G_Domain.map_iZne2cgID[b] = b + 1;

			// Distribute zones evenly trough procs
			const int nAvg = G_Domain.nZones / G_State.mpiNProcs;
			const int nResdl = G_Domain.nZones % G_State.mpiNProcs;
			for (int r = 0; r < G_State.mpiNProcs; ++r)
			{
				const int bs = r * nAvg + ((r<nResdl) ? r : nResdl);
				const int be = bs + nAvg + ((r<nResdl) ? 1 : 0) - 1;

				for (int b = bs; b <= be; ++b)
					G_State.map_zone2rank[b] = r;
			}
		} while (false);

		MPI_Bcast(G_Domain.map_iZne2cgID, G_Domain.nZones, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(G_State.map_zone2rank, G_Domain.nZones, MPI_INT, 0, MPI_COMM_WORLD);

		// Set zone indices range for the current MPI rank
		G_Domain.bs = -1;
		for (int iZne = 0; iZne < G_Domain.nZones; ++iZne)
		{
			const int& r = G_State.map_zone2rank[iZne];

			if (r == G_State.mpiRank)
			{
				if (G_Domain.bs == -1)  G_Domain.bs = iZne;
				G_Domain.be = iZne;
			}
		}

		return true;
}