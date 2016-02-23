//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/parabolic_profile.h"

//libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_function.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exact_solution.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag );

class TurbulentBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  TurbulentBCFactory(libMesh::MeshFunction* _turbulent_bc_values)
    : GRINS::BoundaryConditionsFactory(),
      turbulent_bc_values(_turbulent_bc_values)
  { return; };

  ~TurbulentBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );

private:
  // A pointer to a TurbulentBdyFunction object that build_dirichlet can use to set bcs
  libMesh::MeshFunction* turbulent_bc_values;
};

// Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
class TurbulentBdyFunctionU : public libMesh::FunctionBase<libMesh::Number>
{
public:
  TurbulentBdyFunctionU (libMesh::MeshFunction* _turbulent_bc_values) :
    turbulent_bc_values(_turbulent_bc_values)
  { this->_initialized = true; }

  virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const libMesh::Point& p,
                           const libMesh::Real t,
                           libMesh::DenseVector<libMesh::Number>& output)
  {
    output.resize(1);
    output.zero();

    // Since the turbulent_bc_values object has a solution from a 1-d problem, we have to zero out the y coordinate of p
    libMesh::Point p_copy(p);
    // Also, the 1-d solution provided is on the domain [0, 1] on the x axis and we need to map this to the corresponding point on the y axis
    p_copy(0) = p_copy(1);
    p_copy(1)= 0.0;
    // Also, the 1-d solution provided is actually a symmetry solution, so we have to make the following map
    // x_GRINS < 0.5 => x_meshfunction = 2*x_GRINS , x_GRINS >= 0.5 => x_GRINS = 1 - x_GRINS, x_meshfunction = 2*x_GRINS
    if(p_copy(0) > 0.5)
      {
        p_copy(0) = 1 - p_copy(0);
      }
    p_copy(0) = 2*p_copy(0);

    libMesh::DenseVector<libMesh::Number> u_nu_values;
    turbulent_bc_values->operator()(p_copy, t, u_nu_values);

    output(0) = u_nu_values(0)/21.995539;
  }

  virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
  { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunctionU(turbulent_bc_values)); }

private:
  libMesh::MeshFunction* turbulent_bc_values;
};

// Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
class TurbulentBdyFunctionNu : public libMesh::FunctionBase<libMesh::Number>
{
public:
  TurbulentBdyFunctionNu (libMesh::MeshFunction* _turbulent_bc_values) :
    turbulent_bc_values(_turbulent_bc_values)
  { this->_initialized = true; }

  virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const libMesh::Point& p,
                           const libMesh::Real t,
                           libMesh::DenseVector<libMesh::Number>& output)
  {
    output.resize(1);
    output.zero();

    // Since the turbulent_bc_values object has a solution from a 1-d problem, we have to zero out the y coordinate of p
    libMesh::Point p_copy(p);
    // Also, the 1-d solution provided is on the domain [0, 1] on the x axis and we need to map this to the corresponding point on the y axis
    p_copy(0) = p_copy(1);
    p_copy(1)= 0.0;
    // Also, the 1-d solution provided is actually a symmetry solution, so we have to make the following map
    // x_GRINS < 0.5 => x_meshfunction = 2*x_GRINS , x_GRINS >= 0.5 => x_GRINS = 1 - x_GRINS, x_meshfunction = 2*x_GRINS
    if(p_copy(0) > 0.5)
      {
        p_copy(0) = 1 - p_copy(0);
      }
    p_copy(0) = 2*p_copy(0);

    libMesh::DenseVector<libMesh::Number> u_nu_values;
    turbulent_bc_values->operator()(p_copy, t, u_nu_values);

    output(0) = u_nu_values(1)/(2.0*21.995539);
  }

  virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
  { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunctionNu(turbulent_bc_values)); }

private:
  libMesh::MeshFunction* turbulent_bc_values;
};


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
  libMesh::LibMeshInit libmesh_init(argc, argv);

  // Build a 1-d turbulent_bc_system to get the bc data from files
  libMesh::SerialMesh mesh(libmesh_init.comm());

  GetPot command_line(argc,argv);

  std::string oned_mesh = command_line("mesh-1d", "DIE!");
  mesh.read(oned_mesh);

  // And an EquationSystems to run on it
  libMesh::EquationSystems equation_systems (mesh);

  std::string oned_data = command_line("data-1d", "DIE!");
  equation_systems.read(oned_data, libMesh::XdrMODE::READ,
                        libMesh::EquationSystems::READ_HEADER |
                        libMesh::EquationSystems::READ_DATA |
                        libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  // Get a reference to the system so that we can call update() on it
  libMesh::System & turbulent_bc_system = equation_systems.get_system<libMesh::System>("flow");

  // We need to call update to put system in a consistent state
  // with the solution that was read in
  turbulent_bc_system.update();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Prepare a global solution and a MeshFunction of the Turbulent system
  libMesh::AutoPtr<libMesh::MeshFunction> turbulent_bc_values;

  libMesh::AutoPtr<libMesh::NumericVector<libMesh::Number> > turbulent_bc_soln = libMesh::NumericVector<libMesh::Number>::build(equation_systems.comm());

  std::vector<libMesh::Number> flow_soln;

  turbulent_bc_system.update_global_solution(flow_soln);

  turbulent_bc_soln->init(turbulent_bc_system.solution->size(), true, libMesh::SERIAL);

  (*turbulent_bc_soln) = flow_soln;

  std::vector<unsigned int>turbulent_bc_system_variables;
  turbulent_bc_system_variables.push_back(0);
  turbulent_bc_system_variables.push_back(1);

  turbulent_bc_values = libMesh::AutoPtr<libMesh::MeshFunction>
    (new libMesh::MeshFunction(equation_systems,
			       *turbulent_bc_soln,
			       turbulent_bc_system.get_dof_map(),
			       turbulent_bc_system_variables ));

  turbulent_bc_values->init();

  GRINS::SimulationBuilder sim_builder;

  GRINS::SharedPtr<GRINS::BoundaryConditionsFactory> bc_factory( new TurbulentBCFactory(turbulent_bc_values.get()) );

  sim_builder.attach_bc_factory(bc_factory);

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder,
                           libmesh_init.comm() );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  // Solve
  grins.run();

// Get equation systems to create ExactSolution object
  GRINS::SharedPtr<libMesh::EquationSystems> es = grins.get_equation_system();

  // Create Exact solution object and attach exact solution quantities
  //libMesh::ExactSolution exact_sol(*es);

   // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  libMesh::EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = command_line("soln-data", "DIE!");
  es_ref.read( solution_file, libMesh::XdrMODE::DECODE,
               libMesh::EquationSystems::READ_HEADER |
               libMesh::EquationSystems::READ_DATA |
               libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  exact_sol.attach_reference_solution( &es_ref );

  // Now grab the variables for which we want to compare
  unsigned int n_vars = command_line.vector_variable_size("vars");
  std::vector<std::string> vars(n_vars);
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      vars[v] = command_line("vars", "DIE!", v);
    }

  // Now grab the norms to compute for each variable error
  unsigned int n_norms = command_line.vector_variable_size("norms");
  std::vector<std::string> norms(n_norms);
  for( unsigned int n = 0; n < n_norms; n++ )
    {
      norms[n] = command_line("norms", "DIE!", n);
      if( norms[n] != std::string("L2") &&
          norms[n] != std::string("H1") )
        {
          std::cerr << "ERROR: Invalid norm input " << norms[n] << std::endl
                    << "       Valid values are: L2" << std::endl
                    << "                         H1" << std::endl;
        }
    }

  const std::string& system_name = grins.get_multiphysics_system_name();

  // Now compute error for each variable
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      exact_sol.compute_error(system_name, vars[v]);
    }

  int return_flag = 0;

  double tol = command_line("tol", 1.0e-10);

  // Now test error for each variable, for each norm
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      for( unsigned int n = 0; n < n_norms; n++ )
        {
          test_error_norm( exact_sol, system_name, vars[v], norms[n], tol, return_flag );
        }
    }

  return return_flag;
}

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag )
{
  // We don't set return_flag unless we are setting it 1
  // since this function gets called multiple times and we don't
  // want to overwrite a previous "fail" (return_flag = 1) with
  // a "pass" (return_flag = 0)

  double error = 0.0;

  std::cout << "==========================================================" << std::endl
            << "Checking variable " << var << " using error norm " << norm << " with tol " << tol << "...";

  if( norm == std::string("L2") )
    {
      error = exact_sol.l2_error(system_name, var);
    }
  else if( norm == std::string("H1") )
    {
      error = exact_sol.h1_error(system_name, var);
    }
  else
    {
      std::cerr << "ERROR: Invalid norm " << norm << std::endl;
      exit(1);
    }

  if( error > tol )
    {
      return_flag = 1;

      std::cerr << "Tolerance exceeded for generic regression test!" << std::endl
                << "tolerance     = " << tol << std::endl
                << "norm of error = " << error << std::endl
                << "norm type     = " << norm << std::endl
                << "var           = " << var << std::endl;
    }
  else
    {
      std::cout << "PASSED!" << std::endl
                << "==========================================================" << std::endl;
    }

  return;
}


std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > TurbulentBCFactory::build_dirichlet( )
{
  GRINS::SharedPtr<libMesh::FunctionBase<libMesh::Number> > turbulent_inlet_u( new TurbulentBdyFunctionU(this->turbulent_bc_values) );

  GRINS::SharedPtr<libMesh::FunctionBase<libMesh::Number> > turbulent_inlet_nu( new TurbulentBdyFunctionNu(this->turbulent_bc_values) );

  GRINS::DBCContainer cont_u;
  cont_u.add_var_name( "u" );
  cont_u.add_bc_id( 3 );

  cont_u.set_func( turbulent_inlet_u );

  GRINS::DBCContainer cont_nu;
  cont_nu.add_var_name( "nu" );
  cont_nu.add_bc_id( 3 );

  cont_nu.set_func( turbulent_inlet_nu );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::PhysicsNaming::incompressible_navier_stokes(),  cont_u) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::PhysicsNaming::spalart_allmaras(),  cont_nu) );

  return mymap;
}