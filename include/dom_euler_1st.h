#pragma once

#include "dom_euler_base.h"

class t_DomEu1st : public t_DomEuBase {
public:
	void calcReconstData() {};
	void calcFaceFlux(int iZone, lint iFace);
};
