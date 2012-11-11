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
// $Id: cantera_chem_thermo_test.C 34504 2012-11-10 06:08:51Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <iomanip>
#include "cantera_transport.h"

int main()
{
  GetPot input( "./input_files/cantera_chem_thermo.in" );

  std::vector<std::string> species(5);
  species[0] = input( "Physics/Chemistry/species", "DIE!", 0 );
  species[1] = input( "Physics/Chemistry/species", "DIE!", 1 );
  species[2] = input( "Physics/Chemistry/species", "DIE!", 2 );
  species[3] = input( "Physics/Chemistry/species", "DIE!", 3 );
  species[4] = input( "Physics/Chemistry/species", "DIE!", 4 );

  GRINS::ChemicalMixture chem_mixture(species);

  GRINS::CanteraTransport cantera_trans(input,chem_mixture);

  double T = 1500.0;

  double P = 100000.0;

  std::vector<double> Y(5,0.2);

  GRINS::ReactingFlowCache cache(T,P,Y);

  const double mu = cantera_trans.mu(cache);
  const double k = cantera_trans.k(cache);

  int return_flag = 0;

  const double tol = 1.0e-15;

  std::cout << std::setprecision(16) << std::scientific
	    << "mu = " << mu << std::endl
	    << "k = " << k << std::endl;

  return return_flag;
}
