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
#include "libmesh/mesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exact_solution.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif


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
 
  // Write out this solution to make sure it was read properly
  std::ostringstream file_name_gmv;
  file_name_gmv << "turbulent_bc.gmv";
  
  libMesh::GMVIO(mesh).write_equation_systems
    (file_name_gmv.str(), equation_systems);

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Prepare a global solution and a MeshFunction of the Turbulent system
  libMesh::AutoPtr<libMesh::MeshFunction> turbulent_bc_values;
      
  libMesh::AutoPtr<libMesh::NumericVector<libMesh::Number> > turbulent_bc_soln = libMesh::NumericVector<libMesh::Number>::build(equation_systems.comm());
        
  std::vector<libMesh::Number> flow_soln;

  turbulent_bc_system.update_global_solution(flow_soln);

  std::cout<<"Turbulent system size: "<<turbulent_bc_system.solution->size()<<std::endl;
  
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

  std::tr1::shared_ptr<TurbulentBCFactory> bc_factory( new TurbulentBCFactory(turbulent_bc_values.get()) );

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
  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > turbulent_inlet_u( new TurbulentBdyFunctionU(this->turbulent_bc_values) );

  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > turbulent_inlet_nu( new TurbulentBdyFunctionNu(this->turbulent_bc_values) );

  GRINS::DBCContainer cont_u;
  cont_u.add_var_name( "u" );
  cont_u.add_bc_id( 3 );
    
  cont_u.set_func( turbulent_inlet_u );

  GRINS::DBCContainer cont_nu;
  cont_nu.add_var_name( "nu" );
  cont_nu.add_bc_id( 3 );
    
  cont_nu.set_func( turbulent_inlet_nu );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont_u) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::spalart_allmaras,  cont_nu) );

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

// Class to construct the Dirichlet boundary object and operator for the inlet u velocity and nu profiles
// class TurbulentBdyFunction : public libMesh::FunctionBase<libMesh::Number>
// {
// public:
//   TurbulentBdyFunction (libMesh::MeshFunction* _turbulent_bc_values) :
//     turbulent_bc_values(_turbulent_bc_values)
//   { this->_initialized = true; }

//   virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
//   { libmesh_not_implemented(); }

//   virtual void operator() (const libMesh::Point& p,
//                            const libMesh::Real t,
//                            libMesh::DenseVector<libMesh::Number>& output)
//   {
//     output.resize(4);
//     output.zero();
    
//     // Since the turbulent_bc_values object has a solution from a 1-d problem, we have to zero out the y coordinate of p
//     libMesh::Point p_copy(p);
//     // Also, the 1-d solution provided is on the domain [0, 1] on the x axis and we need to map this to the corresponding point on the y axis
//     p_copy(0) = p_copy(1);
//     p_copy(1)= 0.0;
//     // Also, the 1-d solution provided is actually a symmetry solution, so we have to make the following map
//     // x_GRINS < 0.5 => x_meshfunction = 2*x_GRINS , x_GRINS >= 0.5 => x_GRINS = 1 - x_GRINS, x_meshfunction = 2*x_GRINS
//     if(p_copy(0) > 0.5)
//       {
//        p_copy(0) = 1 - p_copy(0);
//       }
//     p_copy(0) = 2*p_copy(0);
    
//     libMesh::DenseVector<libMesh::Number> u_nu_values;
//     turbulent_bc_values->operator()(p_copy, t, u_nu_values);    
//     //std::cout<<p(1)<<", "<<u_nu_values(0)<<", "<<u_nu_values(1)<<std::endl;    
    
//     output(0) = 1.0; //u_nu_values(0)/21.995539;
//     output(3) = 0.0; //10000*u_nu_values(1)/(2.0*21.995539);
//     }

//   virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
//   { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new TurbulentBdyFunction(turbulent_bc_values)); }

// private:
//   libMesh::MeshFunction* turbulent_bc_values;
// };
  
