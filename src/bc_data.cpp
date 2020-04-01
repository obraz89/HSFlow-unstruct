#include "bc_data.h"

#include "IniFile.hpp"

t_BCList G_BCList;

void t_BCDataInflow::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {}

void t_BCDataOutFlow::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {}

void t_BCDataEulerWall::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {}

void t_BCDataSym::yield(const t_PrimVars& my_pvs, t_PrimVars& opp_pvs) {}