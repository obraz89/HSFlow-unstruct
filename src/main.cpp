#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"

#include "flow_model.h"

#include "flow_params.h"

#include "settings.h"

int main()
{
	load_settings("test_case/main.ini");

	read_cgns_mesh();

	initialize_flow_model();

	G_Domain.makeTimeStep();

	return 0;
}
