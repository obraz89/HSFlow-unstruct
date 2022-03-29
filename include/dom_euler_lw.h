#pragma once

#include "dom_euler_base.h"

class t_DomEuLW : public t_DomEuBase {
public:
	void prepareBeforeTimeMarch() {};
	void calcFaceFlux(int iZone, lint iFace);
};
