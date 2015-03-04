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
//#include "libmesh/exact_solution.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

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

// class ParabolicBCFactory : public GRINS::BoundaryConditionsFactory
// {
// public:

//   ParabolicBCFactory()
//     : GRINS::BoundaryConditionsFactory()
//   { return; };

//   ~ParabolicBCFactory(){return;};

//   std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
// };

// Function to set the Dirichlet boundary function for the inlet u velocity and nu profiles
class BdyFunction : public FunctionBase<Number>
{
public:
  BdyFunction (MeshFunction)
    : MeshFunction.initialize)
  { this->_initialized = true; }

  virtual Number operator() (const Point&, const Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const Point& p,
                           const Real,
                           DenseVector<Number>& output)
  {
    output.resize(2);
    output.zero();
    const Real y=p(1);
    // Set the parabolic inflow boundary conditions at stations 0 & 1
    output(_u_var) = (_sign)*((y-2) * (y-3));
    output(_v_var) = 0;
  }

  virtual AutoPtr<FunctionBase<Number> > clone() const
  { return AutoPtr<FunctionBase<Number> > (new BdyFunction(_u_var, _v_var, _sign)); }

private:
  const unsigned int _u_var, _v_var;
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
  libMesh::SerialMesh mesh(libmesh_init.comm);
  
  libMesh::AutoPtr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  mesh.read("file.xda");
  
  //mesh.all_second_order();
  
  // And an EquationSystems to run on it
  libMesh::EquationSystems equation_systems (mesh);

  libMesh::LinearImplicitSystem & turbulent_bc_system = equation_systems.add_system<libMesh::LinearImplicitSystem>("Turbulent-BC");

  equation_systems.read("sol.xda", READ,
			libMesh::EquationSystems::READ_HEADER |
  			     libMesh::EquationSystems::READ_DATA |
  			     libMesh::EquationSystems::READ_ADDITIONAL_DATA);
 
  // Prepare a global solution and a MeshFunction of the Turbulent system
  libMesh::AutoPtr<MeshFunction> coarse_values;
  libMesh::AutoPtr<libMesh::NumericVector<Number> > turbulent_bc_soln = libMesh::NumericVector<Number>::build(equation_systems.comm());
      
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

  std::tr1::shared_ptr<ParabolicBCFactory> bc_factory( new ParabolicBCFactory );

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

// std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > ParabolicBCFactory::build_dirichlet( )
// {
//   GRINS::DBCContainer cont;
//   cont.add_var_name( "u" );
//   cont.add_bc_id( 1 );
//   cont.add_bc_id( 3 );
  
//   std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > u_func( new GRINS::ParabolicProfile );

//   cont.set_func( u_func );

//   std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

//   mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont) );

//   return mymap;
// }

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
