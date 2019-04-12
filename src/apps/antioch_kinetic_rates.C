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


// GRINS
#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <iomanip>
#include <vector>
#include <fstream>

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_kinetics.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"

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

  GRINS::AntiochMixture<Antioch::CEACurveFit<libMesh::Real> >
    antioch_mixture(input,GRINS::MaterialsParsing::material_name(input,GRINS::PhysicsNaming::reacting_low_mach_navier_stokes()));

  GRINS::AntiochKinetics<Antioch::CEACurveFit<libMesh::Real> > antioch_kinetics( antioch_mixture );

  libMesh::Real T0 = input( "Conditions/T0", 300.0 );
  libMesh::Real T1 = input( "Conditions/T1", 300.0 );
  libMesh::Real T_inc = input( "Conditions/T_increment", 100.0 );


  libMesh::Real p0 = input( "Conditions/pressure", 1.0e5 );

  const unsigned int n_species = antioch_mixture.n_species();

  std::vector<libMesh::Real> Y(n_species);
  if( input.vector_variable_size( "Conditions/mass_fractions" ) != n_species )
    {
      std::cerr << "Error: mass fractions size not consistent with n_species"
                << std::endl;
      libmesh_error();
    }

  for( unsigned int s = 0; s < n_species; s++ )
    {
      Y[s] = input( "Conditions/mass_fractions", 0.0, s );
    }

  libMesh::Real R_mix = antioch_mixture.R_mix(Y);


  std::vector<libMesh::Real> omega_dot(n_species,0.0);

  libMesh::Real T = T0;

  std::ofstream output;
  output.open( "omega_dot.dat", std::ios::trunc );

  output << "# Species names" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      output << antioch_mixture.species_name( s ) << " ";
    }
  output << std::endl;
  output << "# T [K]                  omega_dot [kg/m^3-s]" << std::endl;

  output.close();

  while( T < T1 )
    {
      Antioch::TempCache<libMesh::Real> T_cache(T);

      libMesh::Real rho = p0/(R_mix*T);

      antioch_kinetics.omega_dot( T_cache, rho, Y, omega_dot );

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

#endif //GRINS_HAVE_ANTIOCH
