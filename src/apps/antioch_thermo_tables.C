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
#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/temp_cache.h"
#include "antioch/cea_mixture_ascii_parsing.h"

int main(int argc, char* argv[])
{
  if( argc < 2 )
    {
      std::cerr << "Error: Must specify at least 1 species name!" << std::endl;
      exit(1);
    }

  unsigned n_species = argc-1;
  std::vector<std::string> species_list(n_species);
  for( unsigned int s = 0; s < n_species; s++)
    {
      species_list[s] = std::string(argv[s+1]);
    }

  Antioch::ChemicalMixture<double> chem_mixture( species_list );

  Antioch::CEAThermoMixture<double> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii_default( cea_mixture );

  Antioch::CEAEvaluator<double> cea_thermo( cea_mixture );

  Antioch::StatMechThermodynamics<double> stat_mech_thermo( chem_mixture );

  double T0 = 300.0;
  double delta_T = 100.0;

  std::ofstream output;
  output.open( "thermo.dat", std::ios::trunc );
  output << "# Species names" << std::endl;
  for( unsigned int s = 0; s < n_species; s++)
    {
      output << species_list[s] << " ";
    }
  output << std::endl;
  output << "# T [K]       CEA h_s         StatMech h_s" << std::endl;
  output.close();

  for( unsigned int i = 0; i < 56; i++ )
    {
      double T = T0 + i*delta_T;

      output.open( "thermo.dat", std::ios::app );
      output << std::scientific << std::setprecision(16);
      output << T << " ";

      Antioch::TempCache<double> T_cache(T);

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << cea_thermo.h(T_cache,s) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.h_tot( s, T ) - stat_mech_thermo.h_tot( s, 298.15 )
            + stat_mech_thermo.e_0( s ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.e_tr( s, T ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.e_vib( s, T ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.e_el( s, T ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.e_0( s ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.cv_trans( s ) << " ";
        }

      for( unsigned int s = 0; s < n_species; s++)
        {
          output << std::scientific << std::setprecision(16)
                 << stat_mech_thermo.cv_rot( s ) << " ";
        }

      output << std::endl;
      output.close();
    }
  return 0;
}

#endif // GRINS_HAVE_ANTIOCH
