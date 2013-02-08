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

#include "grins/arrhenius_rate.h"

int main()
{
  const double Cf = 1.4;
  const double eta = 1.2;
  const double Ea = 5.0;

  GRINS::ArrheniusRate arrhenius_rate(Cf,eta,Ea);

  const double T = 1500.1;
  
  const double rate_exact = Cf*std::pow(T,eta)*std::exp(-Ea/T);

  int return_flag = 0;

  double rate = arrhenius_rate(T);

  const double tol = 1.0e-15;

  if( std::fabs( (rate - rate_exact)/rate_exact ) > tol )
    {
      std::cout << "Error: Mismatch in rate values." << std::endl
		<< "rate(T) = " << rate << std::endl
		<< "rate_exact = " << rate_exact << std::endl;

      return_flag = 1;
    }

  std::cout << "Arrhenius rate: " << arrhenius_rate << std::endl;

  return return_flag;
}
