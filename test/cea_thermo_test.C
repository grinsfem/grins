//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// libMesh
#include "getpot.h"

// GRINS
#include "cea_thermo.h"

double cp( double T, double a0, double a1, double a2, 
	   double a3, double a4, double a5, double a6 );

int test_cp( const std::string& species_name, unsigned int species, double cp_exact, double T,
	     const GRINS::CEAThermodynamics& thermo );

int main()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  GRINS::ChemicalMixture chem_mixture( species_str_list );

  GRINS::CEAThermodynamics thermo( GetPot(), chem_mixture );

  const double P = 100000.0;
  const std::vector<double> mass_fractions( 5, 0.2 );
  const double T1 = 190.0;
  const double T2 = 1500.0;
  const double T3 = 10000.0;

  const double R_N2 = GRINS::Constants::R_universal/28.016;
  const double R_exact = GRINS::Constants::R_universal*( 0.2/28.016 + 0.2/32.0 + 0.2/14.008 + 0.2/16.0 + 0.2/30.008 );

  int return_flag = 0;

  // Test N2 cp
  {
    const double cp_N2_1 = R_N2*cp( T1, 2.21037122e+04, -3.81846145e+02, 6.08273815e+00, 
				    -8.53091381e-03,  1.38464610e-05, -9.62579293e-09,  2.51970560e-12);

    const double cp_N2_2 = R_N2*cp( T2, 5.87709908e+05, -2.23924255e+03,  6.06694267e+00,
				    -6.13965296e-04, 1.49179819e-07, -1.92309442e-11, 1.06194871e-15 );

    const double cp_N2_3 = R_N2*cp( T3, 8.30971200e+08, -6.42048187e+05,  2.02020507e+02, -3.06501961e-02,
				    2.48685558e-06, -9.70579208e-11, 1.43751673e-15);
    
    const GRINS::Species species = chem_mixture.species_list()[0];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, 0, cp_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, 0, cp_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, 0, cp_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  return return_flag;
}

int test_cp( const std::string& species_name, unsigned int species, double cp_exact, double T,
	     const GRINS::CEAThermodynamics& thermo )
{
  int return_flag = 0;

  const double tol = 1.0e-15;

  const double cp = thermo.cp(T, species);

  if( std::fabs( (cp_exact - cp)/cp_exact ) > tol )
    {
      std::cerr << "Error: Mismatch in species specific heat." << std::endl
		<< "species = " << species_name << std::endl
		<< "cp       = " << cp << std::endl
		<< "cp_exact = " << cp_exact << std::endl
		<< "T = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

double cp( double T, double a0, double a1, double a2, 
	   double a3, double a4, double a5, double a6 )
{
  if( T < 200.1)
    T = 200.1;

  return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
}
