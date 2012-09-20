//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_unsteady_solver.h"

GRINS::UnsteadySolver::UnsteadySolver( const GetPot& input )
  : Solver(input),
    _theta( input("unsteady-solver/theta", 0.5 ) ),
    _n_timesteps( input("unsteady-solver/n_timesteps", 1 ) ),
    /*! \todo Is this the best default for delta t?*/
    _deltat( input("unsteady-solver/deltat", 0.0 ) )
{
  return;
}

GRINS::UnsteadySolver::~UnsteadySolver()
{
  return;
}

void GRINS::UnsteadySolver::init_time_solver(GRINS::MultiphysicsSystem* system)
{
  libMesh::EulerSolver* time_solver = new libMesh::EulerSolver( *(system) );

  system->time_solver = libMesh::AutoPtr<TimeSolver>(time_solver);

  // Set theta parameter for time-stepping scheme
  time_solver->theta = this->_theta;

  return;
}

void GRINS::UnsteadySolver::solve( GRINS::MultiphysicsSystem* system,
				   std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				   std::tr1::shared_ptr<GRINS::Visualization> vis,
				   bool output_vis,
				   bool output_residual )
{
  system->deltat = this->_deltat;
  
  Real time;

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
    {
      std::cout << "==========================================================" << std::endl
		<< "                 Beginning time step " << t_step  << std::endl
		<< "==========================================================" << std::endl;

      // GRVY timers contained in here (if enabled)
      system->solve();

      time = system->time;

      if( output_vis ) vis->output( equation_system, t_step, time );

      if( output_residual ) vis->output_residual( equation_system, system, t_step, time );

      // Advance to the next timestep
      system->time_solver->advance_timestep();
    }

  return;
}
