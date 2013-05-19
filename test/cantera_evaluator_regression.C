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

//C++
#include <iomanip>

// GRINS
#include "grins_config.h"
#include "grins/cantera_mixture.h"
#include "grins/cantera_evaluator.h"
#include "grins/cached_values.h"

// libMesh
#include "libmesh/getpot.h"

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

  GRINS::CanteraMixture mixture( input );
  GRINS::CanteraEvaluator gas(mixture);

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

  gas.omega_dot( cache, 0, omega_dot );

  const double cv = gas.cv( cache, 0 );
  const double cp = gas.cp( cache, 0 );

  const double mu = gas.mu( cache, 0 );
  const double k = gas.k( cache, 0 );
  
  std::vector<libMesh::Real> D(5,0.0);

  gas.D( cache, 0, D );

  std::vector<double> h(5,0.0);
  gas.h_s( cache, 0, h );

  int return_flag = 0;
  
  double tol = 1.0e-15;
  
  const double cv_reg = 8.8382964243437857e+02;
  if( std::fabs( (cv - cv_reg)/cv_reg ) > tol )
    {
      std::cerr << "Error: Mismatch in internal energy." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "cv = " << cv << std::endl
		<< "cv_reg = " << cv_reg << std::endl;
      return_flag = 1;
    }

  const double cp_reg = 1.2732313697364564e+03;
  if( std::fabs( (cp - cp_reg)/cp_reg ) > tol )
    {
      std::cerr << "Error: Mismatch in internal energy." << std::endl
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

  od_reg[0] = 9.3626353539094969e+04;
  od_reg[1] = -3.3748303628338216e+05;
  od_reg[2] = 2.5813337444763799e+05;
  od_reg[3] = -2.1412192748531760e+05;
  od_reg[4] = 1.9984523578196683e+05;

  for( unsigned int i = 0; i < 5; i++ )
    {
      if( std::fabs( (omega_dot[i] - od_reg[i])/od_reg[i] ) > tol )
	{
	  std::cerr << "Error: Mismatch in internal energy." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "i = " << i << std::endl
		    << "omega_dot = " << omega_dot[i] << std::endl
		    << "od_reg = " << od_reg[i] << std::endl;
	  return_flag = 1;
	}
    }

  const double mu_reg = 5.4816457619629627e-05;
  const double k_reg = 9.2830809315384427e-02;
  std::vector<libMesh::Real> D_reg(5,0.0);
  D_reg[0] = 3.4493708884909694e-04;
  D_reg[1] = 3.3630005706895442e-04;
  D_reg[2] = 3.4039101336823794e-04;
  D_reg[3] = 5.2738514208048540e-04;
  D_reg[4] = 5.1755853748986337e-04;

  if( std::fabs( (mu_reg - mu)/mu ) > tol )
    {
      std::cerr << "Error: Mismatch in viscosity." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "mu     = " << mu << std::endl
		<< "mu_reg = " << mu_reg << std::endl;
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

  std::vector<double> h_reg(5,0.0);
  h_reg[0] = 1.3709248272267890e+06;
  h_reg[1] = 1.2692054328083945e+06;
  h_reg[2] = 4.3659730250572553e+06;
  h_reg[3] = 3.5529883128718123e+07;
  h_reg[4] = 1.7154994250083648e+07;

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

  tol = 6.0e-15;

  if( std::fabs( (k_reg - k)/k ) > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "k     = " << k << std::endl
		<< "k_reg = " << k_reg << std::endl;
      return_flag = 1;
    }
#else // GRINS_HAVE_CANTERA
  // automake expects 77 for a skipped test
  int return_flag = 77;
#endif

  return return_flag;
}
