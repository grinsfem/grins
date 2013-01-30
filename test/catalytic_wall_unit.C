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

// C++
#include <iomanip>

// GRINS
#include "grins/chemical_mixture.h"
#include "grins/catalytic_wall.h"

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

  const unsigned int N_index = 2;

  const GRINS::VariableIndex T_var_dummy = 5;

  const double gamma = 0.03;

  GRINS::CatalyticWall wall_N( chem_mixture, N_index, T_var_dummy, gamma );
  GRINS::CatalyticWall wall_N2( chem_mixture, N_index, T_var_dummy, -gamma );

  const double w_s = 0.2;

  const double rho = 1.0e-3;
  const double rho_s = rho*w_s;

  const double T = 620.1;
  const double R_N = chem_mixture.R( chem_mixture.active_species_name_map().find("N")->second );
  const double M_N = chem_mixture.M( chem_mixture.active_species_name_map().find("N")->second );
  const double R = 30.1;

  const double omega_dot_exact = rho_s*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi*M_N) );
  const double domega_dot_dT_exact = -0.5*rho_s*gamma*std::sqrt( R_N/(T*GRINS::Constants::two_pi*M_N) );
  const double drho_dws = -rho*rho_s/R;
  const double domega_dot_dws_exact = drho_dws*w_s*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi*M_N) )
                                    + rho*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi*M_N) );

  int return_flag = 0;

  const double omega_dot_N = wall_N.omega_dot( rho_s, T );
  const double omega_dot_N2 = wall_N2.omega_dot( rho_s, T );

  const double domega_dot_dT_N = wall_N.domega_dot_dT( rho_s, T );
  const double domega_dot_dT_N2 = wall_N2.domega_dot_dT( rho_s, T );

  const double domega_dot_dws_N = wall_N.domega_dot_dws( rho_s, w_s, T, R );
  const double domega_dot_dws_N2 = wall_N2.domega_dot_dws( rho_s, w_s, T, R );

  const double tol = 1.0e-15;

  /* omega_dot tests */
  {
    double rel_error = std::fabs( (omega_dot_N - omega_dot_exact)/omega_dot_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in omega_dot_N!" << std::endl
		  << "omega_dot_N = " << omega_dot_N << std::endl
		  << "omega_dot_exact = " << omega_dot_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }

  {
    double rel_error = std::fabs( (omega_dot_N2 + omega_dot_exact)/omega_dot_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in omega_dot_N2!" << std::endl
		  << "omega_dot_N2    = " << omega_dot_N2 << std::endl
		  << "omega_dot_exact = " << omega_dot_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }

  /* domega_dot_dT tests */
  {
    double rel_error = std::fabs( (domega_dot_dT_N - domega_dot_dT_exact)/domega_dot_dT_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in domega_dot_dT_N!" << std::endl
		  << "domega_dot_dT_N = " << domega_dot_dT_N << std::endl
		  << "domega_dot_dT_exact = " << domega_dot_dT_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }

  {
    double rel_error = std::fabs( (domega_dot_dT_N2 + domega_dot_dT_exact)/domega_dot_dT_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in domega_dot_dT_N2!" << std::endl
		  << "domega_dot_dT_N2    = " << domega_dot_dT_N2 << std::endl
		  << "domega_dot_dT_exact = " << domega_dot_dT_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }

  /* domega_dot_dws tests */
  {
    double rel_error = std::fabs( (domega_dot_dws_N - domega_dot_dws_exact)/domega_dot_dws_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in domega_dot_dws_N!" << std::endl
		  << "domega_dot_dws_N = " << domega_dot_dws_N << std::endl
		  << "domega_dot_dws_exact = " << domega_dot_dws_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }

  {
    double rel_error = std::fabs( (domega_dot_dws_N2 + domega_dot_dws_exact)/domega_dot_dws_exact );

    if( rel_error > tol )
      {
	std::cerr << "Mismatch in domega_dot_dws_N2!" << std::endl
		  << "domega_dot_dws_N2    = " << domega_dot_dws_N2 << std::endl
		  << "domega_dot_dws_exact = " << domega_dot_dws_exact << std::endl
		  << "rel error = " << rel_error << std::endl;

	return_flag = 1;
      }
  }
    
  return return_flag;
}

