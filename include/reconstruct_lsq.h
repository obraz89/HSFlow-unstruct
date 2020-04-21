#pragma once

#include "mesh.h"

struct t_ReconstDataLSQ {

	// scale used to make reconstruction matrix O(1)
	double r;
	// inverse reconst matrix multiplied by r*r
	t_SqMat3 MInvRR;

};

t_ReconstDataLSQ calcReconstDataLSQ(const t_Cell& cell);
