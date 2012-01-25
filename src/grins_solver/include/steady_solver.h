//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_STEADY_SOLVER_H
#define GRINS_STEADY_SOLVER_H

//GRINS
#include "grins_solver.h"

namespace GRINS
{
  class SteadySolver : public GRINS::Solver
  {
  public:

    SteadySolver();
    virtual ~SteadySolver();

    virtual void solve();
    virtual void init_time_solver();

  };
} // namespace GRINS
#endif // GRINS_STEADY_SOLVER_H
