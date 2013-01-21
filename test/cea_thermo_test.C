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
#include "libmesh/getpot.h"

// GRINS
#include "grins/cea_thermo.h"

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
  const double R_O2 = GRINS::Constants::R_universal/32.0;
  const double R_N = GRINS::Constants::R_universal/14.008;
  const double R_O = GRINS::Constants::R_universal/16.0;
  const double R_NO = GRINS::Constants::R_universal/30.008;
  const double R_exact = mass_fractions[0]*R_N2 + mass_fractions[1]*R_O2 + mass_fractions[2]*R_N +
    mass_fractions[3]*R_O + mass_fractions[4]*R_NO;

  int return_flag = 0;

  // Test N2 cp
  {
    unsigned int index = 0;
    const double cp_N2_1 = R_N2*cp( T1, 2.21037122e+04, -3.81846145e+02, 6.08273815e+00, 
				    -8.53091381e-03,  1.38464610e-05, -9.62579293e-09,  2.51970560e-12);

    const double cp_N2_2 = R_N2*cp( T2, 5.87709908e+05, -2.23924255e+03,  6.06694267e+00,
				    -6.13965296e-04, 1.49179819e-07, -1.92309442e-11, 1.06194871e-15 );

    const double cp_N2_3 = R_N2*cp( T3, 8.30971200e+08, -6.42048187e+05,  2.02020507e+02, -3.06501961e-02,
				    2.48685558e-06, -9.70579208e-11, 1.43751673e-15);
    
    const GRINS::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test O2 cp
  {
    unsigned int index = 1;
    const double cp_1 = R_O2*cp( T1, -3.42556269e+04, 4.84699986e+02, 1.11901159e+00, 
				 4.29388743e-03, -6.83627313e-07, -2.02337478e-09 , 1.03904064e-12 );

    const double cp_2 = R_O2*cp( T2, -1.03793994e+06, 2.34483275e+03, 1.81972949e+00, 1.26784887e-03, 
				 -2.18807142e-07, 2.05372411e-11, -8.19349062e-16 );

    const double cp_3 = R_O2*cp( T3, 4.97515261e+08, -2.86602339e+05, 6.69015464e+01, -6.16971869e-03,  
				 3.01623757e-07, -7.42087888e-12, 7.27744063e-17);
    
    const GRINS::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test N cp
  {
    unsigned int index = 2;
    const double cp_1 = R_N*cp( T1, 0.00000000e+00, 0.00000000e+00, 2.50000000e+00, 0.00000000e+00,
				0.00000000e+00, 0.00000000e+00, 0.00000000e+00);

    const double cp_2 = R_N*cp( T2, 8.87650138e+04, -1.07123150e+02, 2.36218829e+00, 2.91672008e-04,
				-1.72951510e-07, 4.01265788e-11, -2.67722757e-15 );

    const double cp_3 = R_N*cp( T3, 5.47518105e+08, -3.10757498e+05, 6.91678274e+01, -6.84798813e-03,
				3.82757240e-07, -1.09836771e-11, 1.27798602e-16);
    
    const GRINS::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }


  // Test O cp
  {
    unsigned int index = 3;
    const double cp_1 = R_O*cp( T1, -7.95361130e+03 , 1.60717779e+02, 1.96622644e+00, 1.01367031e-03,
				-1.11041542e-06, 6.51750750e-10, -1.58477925e-13 );

    const double cp_2 = R_O*cp( T2, 2.61902026e+05, -7.29872203e+02, 3.31717727e+00, -4.28133436e-04, 
				1.03610459e-07, -9.43830433e-12, 2.72503830e-16 );

    const double cp_3 = R_O*cp( T3, 1.77900426e+08, -1.08232826e+05, 2.81077837e+01, -2.97523226e-03,
				1.85499753e-07, -5.79623154e-12, 7.19172016e-17 );
    
    const GRINS::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }


  // Test NO cp
  {
    unsigned int index = 4;
    const double cp_1 = R_NO*cp( T1, -1.14391658e+04, 1.53646774e+02, 3.43146865e+00, -2.66859213e-03,
				 8.48139877e-06, -7.68511079e-09, 2.38679758e-12 );

    const double cp_2 = R_NO*cp( T2, 2.23903708e+05, -1.28965624e+03, 5.43394039e+00, -3.65605546e-04, 
				 9.88101763e-08, -1.41608327e-11, 9.38021642e-16 );

    const double cp_3 = R_NO*cp( T3, -9.57530764e+08, 5.91243671e+05, -1.38456733e+02, 1.69433998e-02, 
				 -1.00735146e-06, 2.91258526e-11, -3.29511091e-16 );
    
    const GRINS::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag = test_cp( species_name, index, cp_3, T3, thermo );
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
