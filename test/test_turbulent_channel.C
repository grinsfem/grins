//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "libmesh/linear_implicit_system.h"
//#include "libmesh/exact_solution.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif


class TurbulentBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  TurbulentBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~TurbulentBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );

  // A pointer to a TurbulentBdyFunction object that build_dirichlet can use to set bcs
  
};

// Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
class TurbulentBdyFunction : public libMesh::FunctionBase<libMesh::Number>
{
public:
  TurbulentBdyFunction (libMesh::MeshFunction* _turbulent_bc_values) :
    turbulent_bc_values(_turbulent_bc_values)
  { this->_initialized = true; }

  virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const libMesh::Point& p,
                           const libMesh::Real t,
                           libMesh::DenseVector<libMesh::Number>& output)
  {
    output.resize(4);
    output.zero();
    // Since the turbulent_bc_values object has a solution from a 1-d problem, we have to zero out the y coordinate of p
    p(1) = 0.0;
    turbulent_bc_values->operator()(p, t, output);    
  }

  virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
  { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunction(turbulent_bc_values)); }

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
    
  mesh.read("/home/vikram/grins/test/test_data/turbulent_channel_Re944_grid.xda");
  
  //mesh.all_second_order();
  
  // And an EquationSystems to run on it
  libMesh::EquationSystems equation_systems (mesh);

  libMesh::LinearImplicitSystem & turbulent_bc_system = equation_systems.add_system<libMesh::LinearImplicitSystem>("Turbulent-BC");

  equation_systems.read("/home/vikram/grins/test/test_data/turbulent_channel_soln.xda", libMesh::XdrMODE::READ,
			libMesh::EquationSystems::READ_HEADER |
  			     libMesh::EquationSystems::READ_DATA |
  			     libMesh::EquationSystems::READ_ADDITIONAL_DATA);
 
  // Prepare a global solution and a MeshFunction of the Turbulent system
  libMesh::AutoPtr<libMesh::MeshFunction> turbulent_bc_values;
      
libMesh::AutoPtr<libMesh::NumericVector<libMesh::Number> > turbulent_bc_soln = libMesh::NumericVector<libMesh::Number>::build(turbulent_bc_system.comm());
      

  std::vector<unsigned int>turbulent_bc_system_variables;
  turbulent_bc_system_variables.push_back(0);
  turbulent_bc_system_variables.push_back(1);
  
  turbulent_bc_values = libMesh::AutoPtr<libMesh::MeshFunction>
    (new libMesh::MeshFunction(equation_systems,
			       *turbulent_bc_soln,
			       turbulent_bc_system.get_dof_map(),
			       turbulent_bc_system_variables ));
  
  turbulent_bc_values->init();    

  TurbulentBdyFunction turbulent_inlet(turbulent_bc_values.get());

  const libMesh::boundary_id_type left_inlet_id = 0;
  std::set<libMesh::boundary_id_type> left_inlet_bdy;
  left_inlet_bdy.insert(left_inlet_id);

  // The uv identifier for the setting the inlet and wall velocity boundary conditions
  std::vector<unsigned int> unu(1, 0);
  unu.push_back(3);
  
//some_system.get_dof_map().add_dirichlet_boundary
//(libMesh::DirichletBoundary (left_inlet_bdy, unu, &turbulent_inlet));

  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<TurbulentBCFactory> bc_factory( new TurbulentBCFactory );

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
  //std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();

  // Create Exact solution object and attach exact solution quantities
  //libMesh::ExactSolution exact_sol(*es);

  //exact_sol.attach_exact_value(&exact_solution);
  //exact_sol.attach_exact_deriv(&exact_derivative);
  
  // Compute error and get it in various norms
  //exact_sol.compute_error("GRINS", "u");

  //double l2error = exact_sol.l2_error("GRINS", "u");
  //double h1error = exact_sol.h1_error("GRINS", "u");

  // Needs to change to 1 based on comparison
  int return_flag = 0;

  // if( l2error > 1.0e-9 || h1error > 1.0e-9 )
  //   {
  //     return_flag = 1;

  //     std::cout << "Tolerance exceeded for velocity in Poiseuille test." << std::endl
  // 		<< "l2 error = " << l2error << std::endl
  // 		<< "h1 error = " << h1error << std::endl;
  //   }

  // // Compute error and get it in various norms
  // exact_sol.compute_error("GRINS", "p");

  // l2error = exact_sol.l2_error("GRINS", "p");
  // h1error = exact_sol.h1_error("GRINS", "p");

  // if( l2error > 2.0e-9 || h1error > 2.0e-9 )
  //   {
  //     return_flag = 1;

  //     std::cout << "Tolerance exceeded for pressure in Poiseuille test." << std::endl
  // 		<< "l2 error = " << l2error << std::endl
  // 		<< "h1 error = " << h1error << std::endl;
  //   }

  return return_flag;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > TurbulentBCFactory::build_dirichlet( )
{
  GRINS::DBCContainer cont;
  cont.add_var_name( "u" );
  cont.add_bc_id( 1 );
  cont.add_bc_id( 3 );
  
  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > u_func( new GRINS::ParabolicProfile );

  cont.set_func( u_func );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont) );

  return mymap;
}


// libMesh::Number
// exact_solution( const libMesh::Point& p,
// 		const libMesh::Parameters&,   // parameters, not needed
// 		const std::string&,  // sys_name, not needed
// 		const std::string&); // unk_name, not needed);

// libMesh::Gradient
// exact_derivative( const libMesh::Point& p,
// 		  const libMesh::Parameters&,   // parameters, not needed
// 		  const std::string&,  // sys_name, not needed
// 		  const std::string&); // unk_name, not needed);

// libMesh::Number
// exact_solution( const libMesh::Point& p,
// 		const libMesh::Parameters& /*params*/,   // parameters, not needed
// 		const std::string& /*sys_name*/,  // sys_name, not needed
// 		const std::string& var )  // unk_name, not needed);
// {
//   const double x = p(0);
//   const double y = p(1);

//   libMesh::Number f = 0;
//   // Hardcoded to velocity in input file.
//   if( var == "u" ) f = 4*y*(1-y);
//   else if( var == "p" ) f = 120.0 + (80.0-120.0)/5.0*x;
//   else libmesh_assert(false);

//   return f;
// }

// libMesh::Gradient
// exact_derivative( const libMesh::Point& p,
// 		  const libMesh::Parameters& /*params*/,   // parameters, not needed
// 		  const std::string& /*sys_name*/,  // sys_name, not needed
// 		  const std::string& var)  // unk_name, not needed);
// {
//   const double y = p(1);

//   libMesh::Gradient g;

//   // Hardcoded to velocity in input file.
//   if( var == "u" )
//     {
//       g(0) = 0.0;
//       g(1) = 4*(1-y) - 4*y;

// #if LIBMESH_DIM > 2
//       g(2) = 0.0;
// #endif
//     }

//   if( var == "p" )
//     {
//       g(0) = (80.0-120.0)/5.0;
//       g(1) = 0.0;

// #if LIBMESH_DIM > 2
//       g(2) = 0.0;
// #endif
//     }
//   return g;
// }
