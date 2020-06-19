#pragma once

#include "ghost_common.h"

#include "ghost_uvwpt.h"

#include "dom_ns_base.h"

class t_GhostMngNS : public t_GhostMng5 {

	t_DomNSBase* _pDomNS;

public:

	virtual t_Dom5& getDom() { return *_pDomNS; };
	virtual const t_Dom5& getDom() const { return *_pDomNS; };

	void setDom(t_DomNSBase& a_dom) {

		t_GhostMngBase::setDom(a_dom);

		_pDomNS = &a_dom;

	}

};

extern t_GhostMngNS G_GhostMngNS;
