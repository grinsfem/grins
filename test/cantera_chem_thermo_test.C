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
#include "cantera_singleton.h"

int main()
{
  GetPot input( "./input_files/cantera_chem_thermo.in" );

  Cantera::IdealGasMix& cantera = GRINS::CanteraSingleton::cantera_instance( input );

  double T = 1500.0;

  double P = 100000.0;

  std::vector<double> Y(5,0.2);
  //Y[0] = 0.78;
  //Y[1] = 0.22;

  cantera.setState_TPY(T,P,&Y[0]);

  std::vector<double> omega_dot(5,0.0);

  cantera.getNetProductionRates(&omega_dot[0]);
  const double e = cantera.intEnergy_mass();
  const double cv = cantera.cv_mass();
  const double cp = cantera.cp_mass();

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
  od_reg[0] = 3.3421893152544762e+03;
  od_reg[1] = -1.0546740386620191e+04;
  od_reg[2] = 8.6026851320309106e+03;
  od_reg[3] = -1.5287063762539863e+04;
  od_reg[4] = 1.2490795641209472e+04;

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
  

  /*
  std::cout << std::setprecision(16) << std::scientific
	    << "omega_dot = " << omega_dot[0] << ", " << omega_dot[1] << ", " << omega_dot[2]
	    << ", " << omega_dot[3] << ", " << omega_dot[4] << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "e = " << e << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "cv = " << cv << std::endl;
  std::cout << std::setprecision(16) << std::scientific << "cp = " << cp << std::endl;
  */

  return return_flag;
}
