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


// GRINS
#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// C++
#include <fstream>
#include <iomanip>

// GRINS
#include "grins/cantera_mixture.h"
#include "grins/cantera_kinetics.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh_common.h"

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );

  GRINS::CanteraMixture mixture( input );

  GRINS::CanteraKinetics kinetics( mixture );

  libMesh::Real T0 = input( "Conditions/T0", 300.0 );
  libMesh::Real T1 = input( "Conditions/T1", 300.0 );
  libMesh::Real T_inc = input( "Conditions/T_increment", 100.0 );
  
  libMesh::Real rho = input( "Conditions/density", 1.0e-3 );

  const unsigned int n_species = mixture.n_species();

  std::vector<double> Y(n_species,0.0);
  for( unsigned int s = 0; s < n_species; s++ )
    {
      Y[s] = input( "Conditions/mass_fractions", 0.0, s );
    }

  libMesh::Real R_mix = mixture.R_mix(Y);

  std::vector<double> omega_dot(n_species,0.0);

  libMesh::Real T = T0;

  std::ofstream output;
  output.open( "omega_dot.dat", std::ios::trunc );
  
  output << "# Species names" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      output << mixture.species_name( s ) << " ";
    }
  output << std::endl;
  output << "# T [K]                  omega_dot [kg/m^3-s]" << std::endl;

  output.close();

  while( T < T1 )
    { 
      libMesh::Real p = rho*R_mix*T;

      kinetics.omega_dot_TPY( T, p, Y, omega_dot );

      output.open( "omega_dot.dat", std::ios::app );
      output << T << " ";

      for( unsigned int i = 0; i < n_species; i++ )
        {
          output << std::scientific << std::setprecision(16) 
                 << omega_dot[i] << " ";
        }

      output << std::endl;

      output.close();

      T += T_inc;
    }


  return 0;
}

#endif //GRINS_HAVE_CANTERA
