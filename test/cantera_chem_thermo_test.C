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

#include <iomanip>

#include "grins_config.h"
#include "grins/cantera_singleton.h"
#include "grins/cantera_thermo.h"
#include "grins/cantera_kinetics.h"

int main(int argc, char* argv[])
{
#ifdef GRINS_HAVE_CANTERA
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

  Cantera::IdealGasMix& cantera = GRINS::CanteraSingleton::cantera_instance( input );

  GRINS::CanteraThermodynamics cantera_thermo(input,chem_mixture);
  
  GRINS::CanteraKinetics cantera_kinetics(input,chem_mixture);

  double T = 1500.0;

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

  std::vector<double> omega_dot(5,0.0);
  
  cantera_kinetics.omega_dot(cache,0,omega_dot);
  
  const double cv = cantera_thermo.cv( cache, 0 );
  const double cp = cantera_thermo.cp( cache, 0 );

  std::vector<double> h(5,0.0);

  cantera_thermo.h(cache,0,h);

  cantera.setState_TPY(T,P,&Y[0]);
  const double e = cantera.intEnergy_mass();

  int return_flag = 0;
  
  const double tol = 1.0e-15;

  const double e_reg = 1.1354093541825727e+07;
  if( std::fabs( (e - e_reg)/e_reg ) > tol )
    {
      std::cerr << "Error: Mismatch in internal energy." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "e = " << e << std::endl
		<< "e_reg = " << e_reg << std::endl;
      return_flag = 1;
    }
  
  const double cv_reg = 8.8382964243437857e+02;
  if( std::fabs( (cv - cv_reg)/cv_reg ) > tol )
    {
      std::cerr << "Error: Mismatch in cv." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "cv = " << cv << std::endl
		<< "cv_reg = " << cv_reg << std::endl;
      return_flag = 1;
    }

  const double cp_reg = 1.2732313697364564e+03;
  if( std::fabs( (cp - cp_reg)/cp_reg ) > tol )
    {
      std::cerr << "Error: Mismatch in cp." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "cp = " << cp << std::endl
		<< "cp_reg = " << cp_reg << std::endl;
      return_flag = 1;
    }

  std::vector<double> od_reg(5,0.0);
  /* Values before rescaling omega dot by molar mass
  od_reg[0] = 3.3421893152544762e+03;
  od_reg[1] = -1.0546740386620191e+04;
  od_reg[2] = 8.6026851320309106e+03;
  od_reg[3] = -1.5287063762539863e+04;
  od_reg[4] = 1.2490795641209472e+04; */

  od_reg[0] = 9.3634775856169406e+04;
  od_reg[1] = -3.3749569237184612e+05;
  od_reg[2] = 2.5814937544198355e+05;
  od_reg[3] = -2.1414118918565838e+05;
  od_reg[4] = 1.9985273025935155e+05;

  for( unsigned int i = 0; i < 5; i++ )
    {
      if( std::fabs( (omega_dot[i] - od_reg[i])/od_reg[i] ) > tol )
	{
	  std::cerr << "Error: Mismatch in omega_dot." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "i = " << i << std::endl
		    << "omega_dot = " << omega_dot[i] << std::endl
		    << "od_reg = " << od_reg[i] << std::endl;
	  return_flag = 1;
	}
    }
  
  std::vector<double> h_reg(5,0.0);
  h_reg[0] = 1.3708031466651920e+06;
  h_reg[1] = 1.2691593487863187e+06;
  h_reg[2] = 4.3657076051206365e+06;
  h_reg[3] = 3.5526729566942364e+07;
  h_reg[4] = 1.7154371363422986e+07;

  for( unsigned int i = 0; i < 5; i++ )
    {
      if( std::fabs( (h[i] - h_reg[i])/h_reg[i] ) > tol )
	{
	  std::cerr << "Error: Mismatch in internal energy." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "i = " << i << std::endl
		    << "h = " << h[i] << std::endl
		    << "h_reg = " << h_reg[i] << std::endl;
	  return_flag = 1;
	}
    }
  

  /*
  std::cout << std::setprecision(16) << std::scientific
	    << "omega_dot = " << omega_dot[0] << ", " << omega_dot[1] << ", " << omega_dot[2]
	    << ", " << omega_dot[3] << ", " << omega_dot[4] << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "e = " << e << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "cv = " << cv << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "cp = " << cp << std::endl;
  std::cout << std::setprecision(16) << std::scientific
	    << "h = " << h[0] << ", " << h[1] << ", " << h[2]
	    << ", " << h[3] << ", " << h[4] << std::endl;
  */

#else //GRINS_HAVE_CANTERA
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif
  return return_flag;
}
