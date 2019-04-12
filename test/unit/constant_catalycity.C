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
#include "grins/constant_catalycity.h"

int main()
{
  double gamma0 = 0.01;

  GRINS::ConstantCatalycity gamma( gamma0 );

  double T1 = 10.0;

  double tol = std::numeric_limits<double>::epsilon()*10;

  int return_flag = 0;

  if( std::fabs( gamma(T1) - gamma0 ) > tol )
    {
      std::cerr << "Error: mismatch in gamma" << std::endl
                << "       gamma       = " << gamma(T1) << std::endl
                << "       gamma exact = " << gamma0 << std::endl;

      return_flag = 1;
    }

  if( std::fabs( gamma.dT(T1) ) > tol )
    {
      std::cerr << "Error: mismatch in dgamma_dT" << std::endl
                << "       dgamma_dT       = " << gamma.dT(T1) << std::endl
                << "       dgamma_dT exact = " << 0.0 << std::endl;

      return_flag = 1;
    }

  std::vector<double> gamma1;
  gamma1.push_back( 1.0e-5 );

  gamma.set_params( gamma1 );

  if( std::fabs( gamma(T1) - gamma1[0] ) > tol )
    {
      std::cerr << "Error: mismatch in gamma" << std::endl
                << "       gamma       = " << gamma(T1) << std::endl
                << "       gamma exact = " << gamma0 << std::endl;

      return_flag = 1;
    }

  if( std::fabs( gamma.dT(T1) ) > tol )
    {
      std::cerr << "Error: mismatch in dgamma_dT" << std::endl
                << "       dgamma_dT       = " << gamma.dT(T1) << std::endl
                << "       dgamma_dT exact = " << 0.0 << std::endl;

      return_flag = 1;
    }

  return return_flag;
}
