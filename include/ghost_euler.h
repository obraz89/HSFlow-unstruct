#pragma once

#include "ghost_common.h"

#include "dom_euler.h"

class t_GhostMngEuler : public t_GhostMngBase {

	t_DomainEuler* _pDomEu;

public:

	void setDom(t_DomainEuler& a_dom) {

		t_GhostMngBase::setDom(a_dom);

		_pDomEu = &a_dom;

	}

	void exchangeCSV();

};

extern t_GhostMngEuler G_GhostMngEu;
