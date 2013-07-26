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

// This class
#include "grins/mesh_adaptive_solver.h"

// GRINS
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/error_vector.h"
#include "libmesh/steady_solver.h"

GRINS::MeshAdaptiveSolver::MeshAdaptiveSolver( const GetPot& input )
  : Solver( input )
{
  std::string param_path = "Adaptivity/"; // general section for model- and hp-adaptivity
  this->_max_r_steps = input( param_path + "max_r_steps", 5 );
  this->_coarsen_by_parents = true;
  this->_absolute_global_tolerance = input( param_path + "absolute_global_tolerance", 0 );
  this->_nelem_target = input( param_path + "nelem_target", 0 );
  this->_refine_fraction = input( param_path + "refine_percentage", 0.2 );
  this->_coarsen_fraction = input( param_path + "coarsen_percentage", 0.2 );
  this->_coarsen_threshold = input( param_path + "coarsen_threshold", 0 );
  this->_output_adjoint_sol = input( param_path + "output_adjoint_sol", false );

  // Output options
  this->_plot_cell_errors = input( param_path + "plot_cell_errors", false );
  this->_error_plot_prefix = input( param_path + "error_plot_prefix", "cell_error" );

  return;
}

GRINS::MeshAdaptiveSolver::~MeshAdaptiveSolver()
{
  return;
}

void GRINS::MeshAdaptiveSolver::init_time_solver( GRINS::MultiphysicsSystem* system )
{
  libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(system) );

  system->time_solver = libMesh::AutoPtr<libMesh::TimeSolver>( time_solver );

  return;
}

void GRINS::MeshAdaptiveSolver::solve( SolverContext& context )
{
  // Mesh and mesh refinement
  libMesh::MeshBase& mesh = context.equation_system->get_mesh();
  build_mesh_refinement( mesh );

  // This output cannot be toggled in the input file.
  out << "Performing " << this->_max_r_steps << " adaptive refinements" << std::endl;

  // GRVY timers contained in here (if enabled)
  for ( unsigned int r_step = 0;
        r_step < this->_max_r_steps;
        r_step++ )
    {
      // We can't adapt to both a tolerance and a target mesh size
      if( this->_absolute_global_tolerance != 0. )
	{
	  libmesh_assert_equal_to( this->_nelem_target, 0 );
	}

      // If we aren't adapting to a tolerance we need a target mesh size
      else
	{
	  libmesh_assert_greater( _mesh_refinement->nelem_target(), 0 );
	}

      // Solve the forward problem
      context.system->solve();
      libMesh::NumericVector<Number> &primal_solution = *(context.system->solution);
      if( context.output_vis && ! this->_output_adjoint_sol ) context.vis->output( context.equation_system );
    
	  // Make sure we get the contributions to the adjoint RHS from the sides
	  system->assemble_qoi_sides = true;

      // Solve adjoint system: takes transpose of stiffness matrix and solves resulting system
      context.system->adjoint_solve();
      libMesh::NumericVector<Number> &dual_solution = context.system->get_adjoint_solution(0);

    // At the moment output data is overwritten every mesh refinement step
    if( output_vis && this->_output_adjoint_sol )
    {
      // Swap primal and dual to write out dual solution
	    primal_solution.swap( dual_solution );	    
      vis->output( equation_system );
	    primal_solution.swap( dual_solution );	    
    }

      if( context.output_residual )
	{
	  context.vis->output_residual( context.equation_system, context.system );
	}

      // Now we construct the data structures for the mesh refinement process	
      libMesh::ErrorVector error;
	  
      // ``error_estimator'' should have been built in simulation_builder.C
      // and attached to the solver in grins_solver already.
      context.error_estimator->estimate_error( *context.system, error );

      // Plot error vector
      if( this->_plot_cell_errors )
	{
	  error.plot_error( this->_error_plot_prefix+".exo", mesh );
	}
    
	  if( this->_absolute_global_tolerance >= 0. && this->_nelem_target == 0.)
	  {
	    _mesh_refinement->flag_elements_by_error_tolerance( error );
	  }
	  // Adaptively refine based on reaching a target number of elements
	  else
	  {
	    if( mesh.n_active_elem() >= this->_nelem_target )
	    {
	      std::cout<<"We reached the target number of elements."<<std::endl <<std::endl;
	      break;
	    }
	    
	    _mesh_refinement->flag_elements_by_nelem_target( error );
	  }
	  _mesh_refinement->refine_and_coarsen_elements();
    
	  // Dont forget to reinit the system after each adaptive refinement!
	  equation_system->reinit();

    // This output cannot be toggled in the input file.
    std::cout << "Refinement step " << r_step+1 << "/" << this->_max_r_steps
        << ": refined mesh to " << mesh.n_active_elem() << " elements."
        << std::endl << std::endl;
  }

  return;
}

void GRINS::MeshAdaptiveSolver::build_mesh_refinement( libMesh::MeshBase &mesh )
{
  this->_mesh_refinement.reset( new libMesh::MeshRefinement( mesh ) );
  this->_mesh_refinement->coarsen_by_parents() = this->_coarsen_by_parents;
  this->_mesh_refinement->absolute_global_tolerance() = this->_absolute_global_tolerance;
  this->_mesh_refinement->nelem_target() = this->_nelem_target;
  this->_mesh_refinement->refine_fraction() = this->_refine_fraction;
  this->_mesh_refinement->coarsen_fraction() = this->_coarsen_fraction;  
  this->_mesh_refinement->coarsen_threshold() = this->_coarsen_threshold;

  return; 
}
