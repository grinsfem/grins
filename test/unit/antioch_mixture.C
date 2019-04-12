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

// C++
#include <iomanip>
#include <limits>
#include <vector>

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"
#include "grins/antioch_mixture_builder_base.h"

// libMesh
#include "libmesh/getpot.h"

#ifdef GRINS_HAVE_ANTIOCH
int main( int argc, char* argv[] )
#else
  int main()
#endif
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

  GRINS::AntiochMixtureBuilderBase builder;
  std::unique_ptr<GRINS::AntiochMixture<Antioch::CEACurveFit<libMesh::Real> > >
    antioch_ptr = builder.build_antioch_mixture<Antioch::CEACurveFit<libMesh::Real> >(input,"TestMaterial");

  const GRINS::AntiochMixture<Antioch::CEACurveFit<libMesh::Real> > & antioch = *antioch_ptr;

  std::vector<double> mass_fractions( 5, 0.2 );

  // 1.0e-3 converts from kg/kmol -> kg/mol
  const double M_N2 = 14.00800*2*1.0e-3;
  const double M_O2 = 16.0000*2*1.0e-3;
  const double M_N = 14.00800*1.0e-3;
  const double M_O = 16.0000*1.0e-3;
  const double M_NO = 30.00800*1.0e-3;

  double R_exact = Antioch::Constants::R_universal<double>()*( mass_fractions[0]/M_N2
                                                               + mass_fractions[1]/M_O2
                                                               + mass_fractions[3]/M_N
                                                               + mass_fractions[4]/M_O
                                                               + mass_fractions[2]/M_NO );

  double M_exact = 1.0/( mass_fractions[0]/M_N2
                         + mass_fractions[1]/M_O2
                         + mass_fractions[3]/M_N
                         + mass_fractions[4]/M_O
                         + mass_fractions[2]/M_NO );

  std::vector<double> X_exact(5, 0.0);
  X_exact[0] = mass_fractions[0]*M_exact/M_N2;
  X_exact[1] = mass_fractions[1]*M_exact/M_O2;
  X_exact[3] = mass_fractions[3]*M_exact/M_N;
  X_exact[4] = mass_fractions[4]*M_exact/M_O;
  X_exact[2] = mass_fractions[2]*M_exact/M_NO;

  int return_flag = 0;

  const double tol = std::numeric_limits<double>::epsilon()*10;

  if( std::fabs( (antioch.R_mix(mass_fractions) - R_exact)/R_exact) > tol )
    {
      std::cerr << "Error: Mismatch in mixture gas constant." << std::endl
                << std::setprecision(16) << std::scientific
                << "R       = " << antioch.R_mix(mass_fractions) << std::endl
                << "R_exact = " << R_exact <<  std::endl;
      return_flag = 1;
    }

  if( std::fabs( (antioch.M_mix(mass_fractions) - M_exact)/M_exact ) > tol )
    {
      std::cerr << "Error: Mismatch in mixture molar mass." << std::endl
                << std::setprecision(16) << std::scientific
                << "M       = " << antioch.M_mix(mass_fractions) << std::endl
                << "M_exact = " << M_exact << std::endl;
      return_flag = 1;
    }

  std::vector<double> X(5);
  antioch.X( antioch.M_mix(mass_fractions), mass_fractions, X );
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

#else //GRINS_HAVE_ANTIOCH
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif

  return return_flag;
}
