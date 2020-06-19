#pragma once

#include "ghost_common.h"

#include "dom_base_uvwpt.h"

class t_GhostMng5 : public t_GhostMngBase {

public:

	virtual t_Dom5& getDom() = 0;
	virtual const t_Dom5& getDom() const= 0;

	virtual void exchangeCSV();

	virtual void exchangeReconstData();

};
