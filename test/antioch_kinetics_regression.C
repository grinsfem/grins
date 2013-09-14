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


#include "grins_config.h"

// C++
#include <iomanip>
#include <limits>
#include <vector>

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_kinetics.h"

// libMesh
#include "libmesh/getpot.h"

int main( int argc, char* argv[] )
{
#ifdef GRINS_HAVE_ANTIOCH
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );
  
  GRINS::AntiochMixture antioch_mixture(input);

  GRINS::AntiochKinetics antioch_kinetics( antioch_mixture );

  libMesh::Real T = 1000;

  Antioch::TempCache<libMesh::Real> T_cache(T);

  libMesh::Real rho = 1.0e-3;

  const unsigned int n_species = 5;

  std::vector<libMesh::Real> Y(n_species,0.2);

  libMesh::Real R_mix = antioch_mixture.R_mix(Y);

  std::vector<libMesh::Real> omega_dot(n_species,0.0);

  antioch_kinetics.omega_dot( T_cache, rho, R_mix, Y, omega_dot );

  for( unsigned int i = 0; i < n_species; i++ )
    {
      std::cout << std::scientific << std::setprecision(16) 
                << "omega_dot(" << i << ") = " << omega_dot[i] << std::endl;
    }

  int return_flag = 0;
  double tol = std::numeric_limits<double>::epsilon()*10;

  // Check that omega_dot sums to 1
  libMesh::Real sum = 0.0;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      sum += omega_dot[s];
    }

  if( std::fabs(sum) > tol )
    {
      std::cout << "Error: Sum of mass sources not equal to zero!" << std::endl
                << "sum = " << sum << std::endl;
      return_flag = 1;
    }

  // Now check against regression values
  std::vector<libMesh::Real> omega_dot_reg(n_species,0.0);
  omega_dot_reg[0] =  4.3563269373325170e-01;
  omega_dot_reg[1] = -3.6798048754802180e+00;
  omega_dot_reg[2] =  2.9971144322942522e+00;
  omega_dot_reg[3] = -1.8347122381073480e+00;
  omega_dot_reg[4] =  2.0817699875600622e+00;

  for( unsigned int s = 0; s < n_species; s++ )
    {
      if( std::fabs( (omega_dot[s] - omega_dot_reg[s])/omega_dot_reg[s] ) > tol )
	{
	  std::cerr << "Error: Mismatch in omega_dot." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "s = " << s << std::endl
		    << "omega_dot = " << omega_dot[s] << std::endl
		    << "omega_dot_reg = " << omega_dot_reg[s] << std::endl;
	  return_flag = 1;
	}
    }
  
#else //GRINS_HAVE_ANTIOCH
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif

  return return_flag;
}
