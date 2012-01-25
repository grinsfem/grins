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

    SteadySolver( const GetPot& input );
    virtual ~SteadySolver();

    virtual void solve( GRINS::Visualization* vis );
    virtual void init_time_solver();

  };
} // namespace GRINS
#endif // GRINS_STEADY_SOLVER_H
