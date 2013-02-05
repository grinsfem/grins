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

// C++
#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h" 

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

int main(int argc, char* argv[])
{
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
      
      Parameters &params = es->parameters;
      Real T_init = libMesh_inputfile("Physics/LowMachNavierStokes/T0", 0.0);
      Real p0_init = libMesh_inputfile("Physics/LowMachNavierStokes/p0", 0.0);

      Real& dummy_T  = params.set<Real>("T_init");
      dummy_T = T_init;

      Real& dummy_p0 = params.set<Real>("p0_init");
      dummy_p0 = p0_init;

      system.project_solution( initial_values, NULL, params );
    }

  grins.run();

  Real qoi = grins.get_qoi(0);

  // Note that this is a *really* coarse mesh. This is just for testing
  // and not even close to the real QoI for this problem.
  const Real exact_qoi = 4.8158910675853654e+00;

  const Real tol = 5.0e-12;

  int return_flag = 0;

  Real rel_error = std::fabs( (qoi-exact_qoi)/exact_qoi );

  if( rel_error > tol )
    {
      return_flag = 1;
      std::cerr << std::setprecision(16)
		<< std::scientific
		<< "Error: QoI value mismatch." << std::endl
		<< "Computed qoi   = " << qoi << std::endl
		<< "Exact qoi      = " << exact_qoi << std::endl
		<< "Relative error = " << rel_error << std::endl;
    }

  return return_flag;
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
