#pragma once

#pragma once

#include "dom_euler_base.h"

class t_DomEu1stImpl : public t_DomEuBase {
public:
	void calcReconstData() {};
	void calcFaceFlux(int iZone, lint iFace);
	virtual void makeTimeStep();

	// test
	void testPetsc();
};

