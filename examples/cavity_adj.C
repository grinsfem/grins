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
// $Id: cavity.C 36675 2013-02-05 18:29:13Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"

// SolutionHistory Includes  ??ARE THESE NEEDED OR INHERITED ELSEWHERE-VARIS 
#include "libmesh/solution_history.h"
#include "libmesh/memory_solution_history.h"

// Might be already needed, not sure, Varis?
#include "libmesh/adjoint_refinement_estimator.h"

// This is where declare the adjoint refined error estimator. This estimator builds an error bound
// for Q(u) - Q(u_h), by solving the adjoint problem on a finer Finite Element space. For more details
// see the description of the Adjoint Refinement Error Estimator in adjoint_refinement_error_estimator.C
AutoPtr<AdjointRefinementEstimator> build_adjoint_refinement_error_estimator(QoISet &qois)
{
  AutoPtr<AdjointRefinementEstimator> error_estimator;

  std::cout<<"Computing the error estimate using the Adjoint Refinement Error Estimator"<<std::endl<<std::endl;

  AdjointRefinementEstimator *adjoint_refinement_estimator = new AdjointRefinementEstimator;

  error_estimator.reset (adjoint_refinement_estimator);

  adjoint_refinement_estimator->qoi_set() = qois;

  // We enrich the FE space for the dual problem by doing 2 uniform h refinements
  // Maybe we should only do one uniform h refinement?  Varis
  adjoint_refinement_estimator->number_h_refinements = 2;

  return error_estimator;
}


// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

// Function for getting initial adjoint variables (all zero in this example)
Real adjoint_initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

int main(int argc, char* argv[])
{
#ifdef GRINS_USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif

  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  
  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);

  // The Memory Solution History object we will set the system SolutionHistory object to
      MemorySolutionHistory cavity_solution_history(system);
      system.time_solver->set_solution_history(cavity_solution_history);
      
      Parameters &params = es->parameters;
      Real T_init = libMesh_inputfile("Physics/LowMachNavierStokes/T0", 0.0);
      Real p0_init = libMesh_inputfile("Physics/LowMachNavierStokes/p0", 0.0);

      Real& dummy_T  = params.set<Real>("T_init");
      dummy_T = T_init;

      Real& dummy_p0 = params.set<Real>("p0_init");
      dummy_p0 = p0_init;

      system.project_solution( initial_values, NULL, params );
    }
  // Add an adjoint vector, this will be computed after the forward
  // time stepping is complete
  //
  // Tell the library not to save adjoint solutions during the forward
  // solve
  //
  // Tell the library not to project this vector, and hence, memory
  // solution history to not save it.
  //
  // Make this vector ghosted so we can localize it to each element
  // later.
  const std::string & adjoint_solution_name = "cavity_adjoint";
  system.add_vector("cavity_adjoint", false, GHOSTED);
#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  //Run the forward simulation, should work fine with GRINS

  grins.run();


  //Now we have to do the adjoint part manually
  // We will declare an error vector for passing to the adjoint refinement error estimator
  ErrorVector QoI_elementwise_error;
  // And an error vector for accumulating the time contributions(For summing, and using for adaptive refinement on the next run)
  ErrorVector QoI_overtime_element_error;

  // Build an adjoint refinement error estimator object
  // WHAT"S THE CORRECT WAY TO ATTACH THE GRINS QOI?  VARIS
  AutoPtr<AdjointRefinementEstimator> adjoint_refinement_error_estimator =
	  build_adjoint_refinement_error_estimator(qois);


  // Now we will solve the backwards in time adjoint problem
  std::cout << std::endl << "Solving the adjoint problem" << std::endl;

  // We need to tell the library that it needs to project the adjoint, so
  // MemorySolutionHistory knows it has to save it

  // Tell the library to project the adjoint vector, and hence, memory solution history to
  // save it
  system.set_vector_preservation(adjoint_solution_name, true);

  std::cout << "Setting adjoint initial conditions Z("<<system.time<<")"<<std::endl;

  // Need to call adjoint_advance_timestep once for the initial condition setup
  std::cout<<"Retrieving solutions at time t="<<system.time<<std::endl;
  system.time_solver->adjoint_advance_timestep();

  // Output the H1 norm of the retrieved solutions (u^i and u^i+1)
  libMesh::out << "|U(" <<system.time + system.deltat<< ")|= " << system.calculate_norm(*system.solution, 0, H1) << std::endl;

  libMesh::out << "|U(" <<system.time<< ")|= " << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) << std::endl;

  // The first thing we have to do is to apply the adjoint initial
  // condition. The user should supply these. Here they are specified
  // in the functions adjoint_initial_value and adjoint_initial_gradient.  
  // FOR TIME_AVERAGED QOFI CAN USE ZERO BC?
  system.project_vector(adjoint_initial_values, NULL, params, system.get_adjoint_solution(0));
  // Make sure adjoint builds RHS for qoi(or we're going to have a very boring adjoint solution for this problem)
  system.assemble_qoi_sides = true;

  // Since we have specified an adjoint solution for the current time (T), set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the adjoint sensitivity method
  

  system.set_adjoint_already_solved(true);

  libMesh::out << "|Z(" <<system.time<< ")|= " << system.calculate_norm(system.get_adjoint_solution(), 0, H1) << std::endl<<std::endl;

  // write_output(equation_systems, param.n_timesteps, "dual");

  // Now that the adjoint initial condition is set, we will start the
  // backwards in time adjoint integration

  // For loop stepping backwards in time
  for (unsigned int t_step=params.initial_timestep;
       t_step != params.initial_timestep + params.n_timesteps; ++t_step)
    {
      //A pretty update message
      std::cout << " Solving adjoint time step " << t_step << ", time = "
		<< system.time << std::endl;

      // The adjoint_advance_timestep
      // function calls the retrieve function of the memory_solution_history
      // class via the memory_solution_history object we declared earlier.
      // The retrieve function sets the system primal vectors to their values
      // at the current timestep
      std::cout<<"Retrieving solutions at time t="<<system.time<<std::endl;
      system.time_solver->adjoint_advance_timestep();

      // Output the H1 norm of the retrieved solution
      libMesh::out << "|U(" <<system.time + system.deltat << ")|= " << system.calculate_norm(*system.solution, 0, H1) << std::endl;

      libMesh::out << "|U(" <<system.time<< ")|= " << system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) << std::endl;

      system.set_adjoint_already_solved(false);

      system.adjoint_solve();

      //Since we are doing error estimates we now  will compute spatial error estimates using the adjoint weighted 
      //error estimates, using the adjoint weighted residual estimator illustrated in example 4.  This  

      // Estimate the error in each element using the Adjoint Refinement estimator
      // does this include -dU/dt,Z?  probably not, check with Roy, VARIS.  
      adjoint_refinement_error_estimator->estimate_error(system, QoI_elementwise_error);      
      // Now we accumulate this, scaled by deltat, into an overall errorvector, to do adaptivity on
      QoI_overtime_element_error+=.5D0*system.deltat*QoI_elementwise_error;

      // Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator 
      //WE CAN IGNORE THIS ONE-VARIS?

      //system.set_adjoint_already_solved(true);

      libMesh::out << "|Z(" <<system.time<< ")|= "<< system.calculate_norm(system.get_adjoint_solution(), 0, H1) << std::endl << std::endl;

      // Get a pointer to the primal solution vector
      NumericVector<Number> &primal_solution = *system.solution;

      // Get a pointer to the solution vector of the adjoint problem for QoI 0
      NumericVector<Number> &dual_solution_0 = system.get_adjoint_solution(0);

      // Swap the primal and dual solutions so we can write out the adjoint solution
      primal_solution.swap(dual_solution_0);

      //      write_output(equation_systems, params.n_timesteps - (t_step + 1), "dual");

      // Swap back
      primal_solution.swap(dual_solution_0);
    }
    // End adjoint timestep loop

  // My bad C++ here probably.  Need to attach QofI_overtime_elementwise_error so we can use get_global?

  QoI_elementwise_error=QoI_overtime_elementwise_error;

  // Output computed error estimate
  std::cout<< "The computed error in QoI is " << std::setprecision(17)
                 << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) << std::endl;

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}

Real initial_values( const Point&, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;

  if( unknown_name == "T" )
    value = params.get<Real>("T_init");

  else if( unknown_name == "p0" )
    value = params.get<Real>("p0_init");

  else
    value = 0.0;

  return value;
}

Real adjoint_initial_values( const Point&, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;
  return value;
}
