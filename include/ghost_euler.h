#pragma once

#include "ghost_common.h"

#include "dom_euler_base.h"

class t_GhostMngEuler : public t_GhostMngBase {

	t_DomEuBase* _pDomEu;

public:

	void setDom(t_DomEuBase& a_dom) {

		t_GhostMngBase::setDom(a_dom);

		_pDomEu = &a_dom;

	}

	void exchangeCSV();

	void exchangeReconstData();

};

extern t_GhostMngEuler G_GhostMngEu;
