#pragma once

#include "ghost_common.h"

#include "ghost_uvwpt.h"

#include "dom_euler_base.h"

class t_GhostMngEuler : public t_GhostMng5 {

	t_DomEuBase* _pDomEu;

public:

	virtual t_Dom5& getDom() { return *_pDomEu; };
	virtual const t_Dom5& getDom() const { return *_pDomEu; };

	void setDom(t_DomEuBase& a_dom) {

		t_GhostMngBase::setDom(a_dom);

		_pDomEu = &a_dom;

	}

};

extern t_GhostMngEuler G_GhostMngEu;
