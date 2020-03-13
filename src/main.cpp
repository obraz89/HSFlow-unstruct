#include "Mesh-CGNS.h"

#include "logging.h"
#include "common_data.h"

void load_settings() {

	G_Domain.nDim = 3;
};

int main()
{
	load_settings();

	read_cgns_mesh();

	return 0;
}
