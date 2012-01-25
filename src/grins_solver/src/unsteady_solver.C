//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "unsteady_solver.h"

GRINS::UnsteadySolver::UnsteadySolver()
{
  return;
}

GRINS::UnsteadySolver::~UnsteadySolver()
{
  return;
}

void GRINS::UnsteadySolver::read_input_options()
{
  // First call parent class to read in those options.
  GRINS::Solver::read_input_options();

  // Now read options specific for this subclass.
  this->_theta       = input("unsteady-solver/theta", 0.5 );
  this->_n_timesteps = input("unsteady-solver/n_timesteps", 1 );

  /*! \todo Is this the best default for delta t?*/
  this->_deltat      = input("unsteady-solver/deltat", 0.0 );
  
  return;
}

void GRINS::UnsteadySolver::init_time_solver()
{
  libMesh::EulerSolver* time_solver = new libMesh::EulerSolver( *(this->_system) );

  this->_system->time_solver = AutoPtr<TimeSolver>(time_solver);

  // Set theta parameter for time-stepping scheme
  time_solver->theta = this->_theta;

  return;
}

void GRINS::UnsteadySolver::solve()
{
  this->_system->deltat = this->_deltat;

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
    {
      // GRVY timers contained in here (if enabled)
      this->_system->solve();

      // Advance to the next timestep in a transient problem
      this->_system->time_solver->advance_timestep();
    }

  return;
}
