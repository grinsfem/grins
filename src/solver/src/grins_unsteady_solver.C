//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/grins_unsteady_solver.h"

// GRINS
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/euler_solver.h"

namespace GRINS
{

  UnsteadySolver::UnsteadySolver( const GetPot& input )
    : Solver(input),
      _theta( input("unsteady-solver/theta", 0.5 ) ),
      _n_timesteps( input("unsteady-solver/n_timesteps", 1 ) ),
      /*! \todo Is this the best default for delta t?*/
      _deltat( input("unsteady-solver/deltat", 0.0 ) )
  {
    return;
  }

  UnsteadySolver::~UnsteadySolver()
  {
    return;
  }

  void UnsteadySolver::init_time_solver(MultiphysicsSystem* system)
  {
    libMesh::EulerSolver* time_solver = new libMesh::EulerSolver( *(system) );

    system->time_solver = libMesh::AutoPtr<TimeSolver>(time_solver);

    // Set theta parameter for time-stepping scheme
    time_solver->theta = this->_theta;

    return;
  }

  void UnsteadySolver::solve( SolverContext& context )
  {
    libmesh_assert( context.system );

    context.system->deltat = this->_deltat;
  
    Real time;

    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
      {
	std::cout << "==========================================================" << std::endl
		  << "                 Beginning time step " << t_step  << std::endl
		  << "==========================================================" << std::endl;

	// GRVY timers contained in here (if enabled)
	context.system->solve();

	time = context.system->time;

	if( context.output_vis )
	  {
	    context.postprocessing->update_quantities( *(context.equation_system) );
	    context.vis->output( context.equation_system, t_step, time );
	  }

	if( context.output_residual ) context.vis->output_residual( context.equation_system, 
								    context.system, t_step, time );

	// Advance to the next timestep
	context.system->time_solver->advance_timestep();
      }

    return;
  }

} // namespace GRINS
