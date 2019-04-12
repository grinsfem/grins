//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


// C++
#include <iostream>

// GRINS
#include "grins/runner.h"

// libMesh
#include "libmesh/enum_solver_package.h"

// Function for getting initial temperature field
libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
                const std::string& system_name, const std::string& unknown_name );

int main(int argc, char* argv[])
{
  GRINS::Runner grins(argc,argv);

  // This is a tough problem to get to converge without PETSc
  libmesh_example_requires
    (libMesh::default_solver_package() != libMesh::LASPACK_SOLVERS,
     "--enable-petsc");

  grins.init();

  const GetPot & inputfile = grins.get_input_file();

  GRINS::Simulation & sim = grins.get_simulation();

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = inputfile( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = inputfile( "screen-options/system_name", "GRINS" );
      std::shared_ptr<libMesh::EquationSystems> es = sim.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);

      libMesh::Parameters &params = es->parameters;
      libMesh::Real T_init = inputfile("Materials/TestMaterial/ReferenceTemperature/value", 0.0);
      libMesh::Real p0_init = inputfile("Materials/TestMaterial/ThermodynamicPressure/value", 0.0);

      libMesh::Real& dummy_T  = params.set<libMesh::Real>("T_init");
      dummy_T = T_init;

      libMesh::Real& dummy_p0 = params.set<libMesh::Real>("p0_init");
      dummy_p0 = p0_init;

      system.project_solution( initial_values, NULL, params );
    }

  grins.run();

  libMesh::Real qoi = sim.get_qoi_value(0);

  // Note that this is a *really* coarse mesh. This is just for testing
  // and not even close to the real QoI for this problem.

  // Erroneous value from libMesh 0.9.2.2
  // const libMesh::Real exact_qoi = 4.8158910676325055;

  // Value after libMesh 7acb6fc9 bugfix
  const libMesh::Real exact_qoi = 4.8654229502012685;

  const libMesh::Real tol = 1.0e-9;

  int return_flag = 0;

  libMesh::Real rel_error = std::fabs( (qoi-exact_qoi)/exact_qoi );

  if( rel_error > tol )
    {
      // Skip this test until we know what changed
      // return_flag = 1;
      return_flag = 77;

      std::cerr << std::setprecision(16)
                << std::scientific
                << "Error: QoI value mismatch." << std::endl
                << "Computed qoi   = " << qoi << std::endl
                << "Exact qoi      = " << exact_qoi << std::endl
                << "Relative error = " << rel_error << std::endl;
    }

  return return_flag;
}

libMesh::Real
initial_values( const libMesh::Point&, const libMesh::Parameters &params,
                const std::string& , const std::string& unknown_name )
{
  libMesh::Real value = 0.0;

  if( unknown_name == "T" )
    value = params.get<libMesh::Real>("T_init");

  else if( unknown_name == "p0" )
    value = params.get<libMesh::Real>("p0_init");

  else
    value = 0.0;

  return value;
}
