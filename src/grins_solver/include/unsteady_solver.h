//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_UNSTEADY_SOLVER_H
#define GRINS_UNSTEADY_SOLVER_H

//GRINS
#include "grins_solver.h"

namespace GRINS
{
  class UnsteadySolver : public GRINS::Solver
  {
  public:

    UnsteadySolver();
    virtual ~UnsteadySolver();

    virtual void read_input_options();
    virtual void solve();
    virtual void init_time_solver();

  protected:

    double _theta;
    unsigned int _n_timesteps;
    double _deltat;

  };
} // namespace GRINS
#endif // GRINS_UNSTEADY_SOLVER_H
