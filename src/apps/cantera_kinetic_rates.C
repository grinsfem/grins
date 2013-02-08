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

// GRINS
#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// Cantera
#include "grins/cantera_kinetics.h"

int main(int argc, char* argv[])
{
  if( argc < 3 )
    {
      std::cerr << "Error: Must specify mixture and input file." << std::endl;
      exit(1);
    }

  std::string mixture = argv[1];
  std::string chem_file = argv[2];

  Cantera::IdealGasMix gas( chem_file, mixture);
  int n_species = gas.nSpecies();

  // 0th argument + input files + P + T + n_species mass_fractions
  if( argc < 3+2+n_species )
    {
      std::cerr << "Error: Expecting mixture, chem_file, P, T, and " << n_species << "mass_fractions" << std::endl;
      exit(1);
    }
  double P = atof(argv[3]);
  double T = atof(argv[4]);

  std::vector<double> Y(n_species,0.0);
  for( int s = 0; s < n_species; s++ )
    {
      Y[s] = atof(argv[5+s]);
    }

  try
    {
      gas.setState_TPY(T,P,&Y[0]);
    }
  catch(Cantera::CanteraError)
    {
      Cantera::showErrors(std::cerr);
      exit(1);
    }
  std::vector<double> omega_dot(n_species,0.0);
  std::vector<double> kfwd(n_species,0.0);
  std::vector<double> kbwd(n_species,0.0);

  try
    {
      gas.getNetProductionRates(&omega_dot[0]);
      gas.getFwdRateConstants(&kfwd[0]);
      gas.getRevRateConstants(&kbwd[0]);
    }
  catch(Cantera::CanteraError)
    {
      Cantera::showErrors(std::cerr);
      exit(1);
    }

  for( int s = 0; s < n_species; s++ )
    {
      std::cout << "omega_dot[" << gas.speciesName(s) << "] = " << omega_dot[s] << std::endl;
    }

  for( unsigned int r = 0; r < gas.nReactions(); r++ )
    {
      std::cout << "kfwd[" << r << "] = " << kfwd[r] << std::endl;
    }

  for( unsigned int r = 0; r < gas.nReactions(); r++ )
    {
      std::cout << "kbwd[" << r << "] = " << kbwd[r] << std::endl;
    }

  return 0;
}

#endif //GRINS_HAVE_CANTERA
