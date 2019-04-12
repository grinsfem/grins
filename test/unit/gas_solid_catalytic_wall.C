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


// C++
#include <iomanip>

// GRINS
#include "grins/math_constants.h"
#include "grins/gas_solid_catalytic_wall.h"
#include "grins/constant_catalycity.h"
#include "grins/cantera_mixture.h"
#include "grins/antioch_chemistry.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"
#include "grins/chemistry_builder.h"

// libMesh
#include "libmesh/getpot.h"

template<typename ChemicalMixture>
int test( ChemicalMixture& chem_mixture )
{
  const unsigned int N_index = chem_mixture.species_index("N");
  const unsigned int C_index = chem_mixture.species_index("C");
  const unsigned int CN_index = chem_mixture.species_index("CN");

  const double gamma = 0.03;

  GRINS::ConstantCatalycity gamma_r( gamma );

  GRINS::GasSolidCatalyticWall<ChemicalMixture> wall( chem_mixture, gamma_r, N_index, C_index, CN_index );

  const double w_s = 0.2;

  const double rho = 1.0e-3;
  const double rho_s = rho*w_s;

  const double T = 620.1;
  const double R_N = chem_mixture.R( chem_mixture.species_index("N") );

  const double R = 30.1;

  const double omega_dot_exact = rho_s*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi) );
  const double domega_dot_dT_exact = 0.5*rho_s*gamma*std::sqrt( R_N/(T*GRINS::Constants::two_pi) );
  const double drho_dws = -rho*rho_s/R;
  const double domega_dot_dws_exact = drho_dws*w_s*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi) )
    + rho*gamma*std::sqrt( R_N*T/(GRINS::Constants::two_pi) );

  const double mdot_C_exact = -omega_dot_exact*chem_mixture.M( chem_mixture.species_index("C") )/chem_mixture.M( chem_mixture.species_index("N") );

  const double mdot_CN_exact = omega_dot_exact*chem_mixture.M( chem_mixture.species_index("CN") )/chem_mixture.M( chem_mixture.species_index("N") );

  int return_flag = 0;

  const double omega_dot_N = wall.omega_dot( rho_s, T );

  const double domega_dot_dT_N = wall.domega_dot_dT( rho_s, T );

  const double domega_dot_dws_N = wall.domega_dot_dws( rho_s, w_s, T, R );

  const double mdot_N = wall.compute_reactant_gas_mass_flux( rho, w_s, T );

  const double mdot_C = wall.compute_reactant_solid_mass_consumption( rho, w_s, T );

  const double mdot_CN = wall.compute_product_mass_flux( rho, w_s, T );

  const double tol = 1.0e-15;

  /* omega_dot tests */
  {
    double rel_error = std::fabs( (omega_dot_N - omega_dot_exact)/omega_dot_exact );

    if( rel_error > tol )
      {
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in omega_dot_N!" << std::endl
                  << "omega_dot_N = " << omega_dot_N << std::endl
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
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in domega_dot_dT_N!" << std::endl
                  << "domega_dot_dT_N = " << domega_dot_dT_N << std::endl
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
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in domega_dot_dws_N!" << std::endl
                  << "domega_dot_dws_N = " << domega_dot_dws_N << std::endl
                  << "domega_dot_dws_exact = " << domega_dot_dws_exact << std::endl
                  << "rel error = " << rel_error << std::endl;

        return_flag = 1;
      }
  }

  /* mdot_N tests */
  {
    double rel_error = std::fabs( (mdot_N-(-omega_dot_exact))/(-omega_dot_exact) );

    if( rel_error > tol )
      {
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in mdot_N!" << std::endl
                  << "mdot_N = " << mdot_N << std::endl
                  << "mdot_N_exact = " << -omega_dot_exact << std::endl
                  << "rel error = " << rel_error << std::endl;

        return_flag = 1;
      }
  }

  /* mdot_C tests */
  {
    double rel_error = std::fabs( (mdot_C - mdot_C_exact)/mdot_C_exact );

    if( rel_error > tol )
      {
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in mdot_C!" << std::endl
                  << "mdot_C = " << mdot_C << std::endl
                  << "mdot_C_exact = " << mdot_C_exact << std::endl
                  << "rel error = " << rel_error << std::endl;

        return_flag = 1;
      }
  }

  /* mdot_CN tests */
  {
    double rel_error = std::fabs( (mdot_CN - mdot_CN_exact)/mdot_CN_exact );

    if( rel_error > tol )
      {
        std::cerr << std::setprecision(16) << std::scientific
                  << "Mismatch in mdot_CN!" << std::endl
                  << "mdot_CN = " << mdot_CN << std::endl
                  << "mdot_CN_exact = " << mdot_CN_exact << std::endl
                  << "rel error = " << rel_error << std::endl;

        return_flag = 1;
      }
  }

  return return_flag;
}


int main(int argc, char* argv[])
{
  if( argc != 3 )
    {
      std::cerr << "Error: must specify the test type (cantera or antioch) and the input file name"
                << std::endl;
      exit(1);
    }

  std::string test_type = argv[1];

  GetPot input( argv[2] );

  int return_flag = 0;

  GRINS::ChemistryBuilder chem_builder;

  if( test_type == "cantera" )
    {
#ifdef GRINS_HAVE_CANTERA
      std::unique_ptr<GRINS::CanteraMixture> chem_uptr;
      chem_builder.build_chemistry(input,"CanteraMaterial",chem_uptr);
      return_flag = test<GRINS::CanteraMixture>( *chem_uptr );
#else
      return_flag = 77;
#endif
    }
  else if( test_type == "antioch" )
    {
#ifdef GRINS_HAVE_ANTIOCH
      std::unique_ptr<GRINS::AntiochChemistry> chem_uptr;
      chem_builder.build_chemistry(input,"AntiochMaterial",chem_uptr);
      return_flag = test<GRINS::AntiochChemistry>( *chem_uptr );
#else
      return_flag = 77;
#endif
    }
  else
    {
      std::cerr << "Invalid test type " << test_type << std::endl;
      return_flag = 1;
    }

  return return_flag;
}
