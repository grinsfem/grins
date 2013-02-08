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

// GRINS
#include "grins/chemical_mixture.h"
#include "grins/cea_thermo.h"
#include "grins/grins_kinetics.h"

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  GetPot input( argv[1] );

  std::vector<std::string> species(5);
  species[0] = input( "Physics/Chemistry/species", "DIE!", 0 );
  species[1] = input( "Physics/Chemistry/species", "DIE!", 1 );
  species[2] = input( "Physics/Chemistry/species", "DIE!", 2 );
  species[3] = input( "Physics/Chemistry/species", "DIE!", 3 );
  species[4] = input( "Physics/Chemistry/species", "DIE!", 4 );

  GRINS::ChemicalMixture chem_mixture(species);

  GRINS::Kinetics kinetics(input,chem_mixture);

  GRINS::CEAThermodynamics thermo(input,chem_mixture);

  const double T = 1500.0;
  const double P = 100000.0;
  const std::vector<double> Y(5,0.2);
  const double R_mix = chem_mixture.R(Y);
  const double rho = P/(R_mix*T);
  
  std::cout << "R_mix = " << R_mix << std::endl;
  std::cout << "rho = " << rho << std::endl;
  
  std::vector<double> molar_densities(5,0.0);
  chem_mixture.molar_densities( rho, Y, molar_densities );

  std::vector<double> h(5,0.0);
  thermo.h_RT_minus_s_R(T,h);

  GRINS::CachedValues cache;

  cache.add_quantity(GRINS::Cache::TEMPERATURE);
  cache.add_quantity(GRINS::Cache::MIXTURE_DENSITY);
  cache.add_quantity(GRINS::Cache::MIXTURE_GAS_CONSTANT);
  cache.add_quantity(GRINS::Cache::MASS_FRACTIONS);
  cache.add_quantity(GRINS::Cache::MOLAR_DENSITIES);
  cache.add_quantity(GRINS::Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY);

  std::vector<double> Tqp(1,T);
  std::vector<double> rhoqp(1,rho);
  std::vector<double> Rqp(1,R_mix);
  std::vector<std::vector<double> > Yqp(1,Y);
  std::vector<std::vector<double> > mqp(1,molar_densities);
  std::vector<std::vector<double> > hqp(1,h);

  cache.set_values(GRINS::Cache::TEMPERATURE, Tqp);
  cache.set_values(GRINS::Cache::MIXTURE_DENSITY, rhoqp);
  cache.set_values(GRINS::Cache::MIXTURE_GAS_CONSTANT, Rqp);
  cache.set_vector_values(GRINS::Cache::MOLAR_DENSITIES, Yqp);
  cache.set_vector_values(GRINS::Cache::MASS_FRACTIONS, mqp);
  cache.set_vector_values(GRINS::Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY, hqp);

  std::vector<double> omega_dot(5,0.0);
  
  kinetics.omega_dot(cache,0,omega_dot);

  std::cout << kinetics.reaction_set() << std::endl;

  int return_flag = 0;

  
  std::cout << std::setprecision(16) << std::scientific
	    << "omega_dot = " << omega_dot[0] << ", " << omega_dot[1] << ", " << omega_dot[2]
	    << ", " << omega_dot[3] << ", " << omega_dot[4] << std::endl;

  return return_flag;
}
