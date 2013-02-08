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

//C++
#include <iomanip>

//GRINS
#include "grins_config.h"
#include "grins/cantera_transport.h"

int main()
{
#ifdef GRINS_HAVE_CANTERA
  GetPot input( "./input_files/cantera_transport.in" );

  std::vector<std::string> species(5);
  species[0] = input( "Physics/Chemistry/species", "DIE!", 0 );
  species[1] = input( "Physics/Chemistry/species", "DIE!", 1 );
  species[2] = input( "Physics/Chemistry/species", "DIE!", 2 );
  species[3] = input( "Physics/Chemistry/species", "DIE!", 3 );
  species[4] = input( "Physics/Chemistry/species", "DIE!", 4 );

  GRINS::ChemicalMixture chem_mixture(species);

  GRINS::CanteraTransport cantera_trans(input,chem_mixture);

  double T = 1000.0;

  double P = 100000.0;

  std::vector<double> Y(5,0.2);

  GRINS::CachedValues cache;

  cache.add_quantity(GRINS::Cache::TEMPERATURE);
  cache.add_quantity(GRINS::Cache::THERMO_PRESSURE);
  cache.add_quantity(GRINS::Cache::MASS_FRACTIONS);

  std::vector<double> Tqp(1,T);
  std::vector<double> Pqp(1,P);
  std::vector<std::vector<double> > Yqp(1,Y);

  cache.set_values(GRINS::Cache::TEMPERATURE, Tqp);
  cache.set_values(GRINS::Cache::THERMO_PRESSURE, Pqp);
  cache.set_vector_values(GRINS::Cache::MASS_FRACTIONS, Yqp);

  const double mu = cantera_trans.mu(cache,0);
  const double k = cantera_trans.k(cache,0);
  
  std::vector<libMesh::Real> D(5,0.0);

  cantera_trans.D(cache, 0, D);

  int return_flag = 0;

  const double tol = 6.0e-15;

  const double mu_reg = 4.2134235819759682e-05;
  const double k_reg = 5.7138665373733508e-02;
  std::vector<libMesh::Real> D_reg(5,0.0);
  D_reg[0] = 1.7611544183904180e-04;
  D_reg[1] = 1.7169898621123060e-04;
  D_reg[2] = 1.7379080956310527e-04;
  D_reg[3] = 2.6991049091576078e-04;
  D_reg[4] = 2.6488528311729070e-04;

  if( std::fabs( (mu_reg - mu)/mu ) > tol )
    {
      std::cerr << "Error: Mismatch in viscosity." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "mu     = " << mu << std::endl
		<< "mu_reg = " << mu_reg << std::endl;
      return_flag = 1;
    }

  if( std::fabs( (k_reg - k)/k ) > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "k     = " << k << std::endl
		<< "k_reg = " << k_reg << std::endl;
      return_flag = 1;
    }

  for( unsigned int i = 0; i < 5; i++ )
    {
      if( std::fabs( (D_reg[i] - D[i])/D[i] ) > tol )
	{
	  std::cerr << "Error: Mismatch in diffusion coefficient." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "i = " << i << std::endl
		    << "D     = " << D[i] << std::endl
		    << "D_reg = " << D_reg[i] << std::endl;
	  return_flag = 1;
	}
    }
  /*
  std::cout << std::setprecision(16) << std::scientific
	    << "mu = " << mu << std::endl
	    << "k = " << k << std::endl;
  for( unsigned int i = 0; i < 5; i++ )
    {
      std::cout << std::setprecision(16) << std::scientific
		<< "D[" << i << "] = " << D[i]
		<< std::endl;
    }
  */
#else //GRINS_HAVE_CANTERA
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif

  return return_flag;
}
