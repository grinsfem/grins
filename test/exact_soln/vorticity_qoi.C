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


#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/runner.h"

//libMesh
#include "libmesh/exact_solution.h"

int main(int argc, char* argv[])
{
  GRINS::Runner grins(argc,argv);
  grins.init();

  // Solve
  grins.run();

  GRINS::Simulation & sim = grins.get_simulation();

  libMesh::Number qoi = sim.get_qoi_value( 0 );

  int return_flag = 0;
  const libMesh::Number exact_value = -0.5;
  const libMesh::Number rel_error = std::fabs( (qoi - exact_value )/exact_value );
  const libMesh::Number tol = 1.0e-11;
  if( rel_error > tol )
    {
      std::cerr << "Computed voriticity QoI mismatch greater than tolerance." << std::endl
                << "Computed value = " << qoi << std::endl
                << "Exact value = " << exact_value << std::endl
                << "Relative error = " << rel_error << std::endl
                << "Tolerance = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}
