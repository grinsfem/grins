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
#include <iomanip>
#include <limits>
#include <vector>

// GRINS
#include "grins/cantera_chemistry.h"

// libMesh
#include "libmesh/getpot.h"

int main( int argc, char* argv[] )
{
#ifdef GRINS_HAVE_CANTERA
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );

  GRINS::CanteraChemistry cantera(input);

  std::vector<double> mass_fractions( 5, 0.2 );

  double R_exact = Cantera::GasConstant*( 0.2/28.016 + 0.2/32.0 + 0.2/14.008 + 0.2/16.0 + 0.2/30.008 );

  double M_exact = 1.0/( 0.2*( 1.0/28.016 + 1.0/32.0 + 1.0/14.008 + 1.0/16.0 + 1.0/30.008) );
  
  std::vector<double> X_exact(5, 0.0);
  X_exact[0] = 0.2*M_exact/28.016;
  X_exact[1] = 0.2*M_exact/32.0;
  X_exact[2] = 0.2*M_exact/14.008;
  X_exact[3] = 0.2*M_exact/16.0;
  X_exact[4] = 0.2*M_exact/30.008;

  int return_flag = 0;

  const double tol = std::numeric_limits<double>::epsilon()*10;

  if( std::fabs( (cantera.R_mix(mass_fractions) - R_exact)/R_exact) > tol )
    {
      std::cerr << "Error: Mismatch in mixture gas constant." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "R       = " << cantera.R_mix(mass_fractions) << std::endl
		<< "R_exact = " << R_exact <<  std::endl;
      return_flag = 1;
    }

  if( std::fabs( (cantera.M_mix(mass_fractions) - M_exact)/M_exact ) > tol )
    {
      std::cerr << "Error: Mismatch in mixture molar mass." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "M       = " << cantera.M_mix(mass_fractions) << std::endl
		<< "M_exact = " << M_exact << std::endl;
      return_flag = 1;
    }
  
  std::vector<double> X(5);
  cantera.X( cantera.M_mix(mass_fractions), mass_fractions, X );
  for( unsigned int s = 0; s < 5; s++ )
    {
      if( std::fabs( (X[s] - X_exact[s])/X_exact[s]) > tol )
	{
	  std::cerr << "Error: Mismatch in mole fraction for species " << s << std::endl
		    << std::setprecision(16) << std::scientific
		    << "X       = " << X[s] << std::endl
		    << "X_exact = " << X_exact[s] << std::endl;
	  return_flag = 1;
	}
    }

#else //GRINS_HAVE_CANTERA
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif

  return return_flag;
}
