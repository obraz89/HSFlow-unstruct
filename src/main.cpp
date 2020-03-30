#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"

#include "flow_model.h"

#include "flow_params.h"

void load_settings() {

	G_Domain.nDim = 3;
};

int main()
{
	load_settings();

	read_cgns_mesh();

	initialize_flow_model();

	initialize_flow_params();

	G_Domain.makeTimeStep();

	return 0;
}
