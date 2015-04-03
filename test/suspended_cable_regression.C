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
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/math_constants.h"

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/exact_solution.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// SuspendedCable
//#include "suspended_cable_solver_factory.h"

// Function for getting initial temperature field
libMesh::Real initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
		                      const std::string& system_name, const std::string& unknown_name );

int main(int argc, char* argv[])
{
#ifdef GRINS_USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif
	// Check command line count.
	if( argc < 3 )
	{
		// TODO: Need more consistent error handling.
		std::cerr << "Error: Must specify libMesh input file." << std::endl;
		exit(1); // TODO: something more sophisticated for parallel runs?
	}

	// libMesh input file should be first argument
	std::string libMesh_input_filename = argv[1];

	// Create our GetPot object.
	GetPot libMesh_inputfile( libMesh_input_filename );

	// GetPot doesn't throw an error for a nonexistent file?
	{
		std::ifstream i(libMesh_input_filename.c_str());
		if (!i)
		{
			std::cerr << "Error: Could not read from libMesh input file "
					<< libMesh_input_filename << std::endl;
			exit(1);
		}
	}

	// Initialize libMesh library.
	libMesh::LibMeshInit libmesh_init(argc, argv);

	libMesh::out << "Starting GRINS with command:\n";
	for (int i=0; i != argc; ++i)
		libMesh::out << argv[i] << ' ';
	libMesh::out << std::endl;

	GRINS::SimulationBuilder sim_builder;

	GRINS::Simulation grins( libMesh_inputfile,
						     sim_builder,
						     libmesh_init.comm() );

	std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );

	// Get equation systems
	std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
	const libMesh::System& system = es->get_system(system_name);

	libMesh::Parameters &params = es->parameters;

	system.project_solution( initial_values, NULL, params );

	grins.run();

	//es->write("suspended_cable_test.xdr");

	// Create Exact solution object and attach exact solution quantities
	libMesh::ExactSolution exact_sol(*es);

	libMesh::EquationSystems es_ref( es->get_mesh() );

	// Filename of file where comparison solution is stashed
	std::string solution_file = std::string(argv[2]);
	es_ref.read( solution_file );

	exact_sol.attach_reference_solution( &es_ref );

	// Compute error and get it in various norms
	exact_sol.compute_error(system_name, "u");
	exact_sol.compute_error(system_name, "v");
	exact_sol.compute_error(system_name, "w");

	double u_l2error = exact_sol.l2_error(system_name, "u");
	double u_h1error = exact_sol.h1_error(system_name, "u");

	double v_l2error = exact_sol.l2_error(system_name, "v");
	double v_h1error = exact_sol.h1_error(system_name, "v");

	double w_l2error = exact_sol.l2_error(system_name, "w");
	double w_h1error = exact_sol.h1_error(system_name, "w");

	int return_flag = 0;

	double tol = 5.0e-8;

	if( u_l2error > tol   || u_h1error > tol   ||
	    v_l2error > tol   || v_h1error > tol   ||
	    w_l2error > tol   || w_h1error > tol     )
	{
	  return_flag = 1;

	  std::cout << "Tolerance exceeded for suspended cable test." << std::endl
		<< "tolerance     = " << tol << std::endl
		<< "u l2 error    = " << u_l2error << std::endl
		<< "u h1 error    = " << u_h1error << std::endl
		<< "v l2 error    = " << v_l2error << std::endl
		<< "v h1 error    = " << v_h1error << std::endl
		<< "w l2 error    = " << w_l2error << std::endl
		<< "w h1 error    = " << w_h1error << std::endl;
	}

	return return_flag;
}

libMesh::Real initial_values( const libMesh::Point& p, const libMesh::Parameters &/*params*/,
		     const std::string& , const std::string& unknown_name )
{
	libMesh::Real value = 0.0;


	if( unknown_name == "u" )
	{
		value = -35*sin(GRINS::Constants::pi*p(0)/400.0);
	}
	else if( unknown_name == "w" )
	{
		value = -55*sin(GRINS::Constants::pi*p(0)/200.0);
	}
	else
	{
		value = 0.0;
	}

	return value;
}
