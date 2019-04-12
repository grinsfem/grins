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
#include "grins/gaussian_xy_profile.h"

// libMesh
#include "libmesh/point.h"

int main()
{
  const double a = 5.0;
  const double mu = 1.5;
  const double sigma = 2.1;
  const double b = 3.2;
  GRINS::GaussianXYProfile profile( a, mu, sigma, b );

  const double x = 1.1;
  const double y = 1.5;

  libMesh::Point p( x, y );

  const double value = profile( p, 0.0 );

  const double r = std::sqrt( x*x + y*y);
  const double exact_value = a*std::exp( -(r-mu)*(r-mu)/(2.0*sigma*sigma) ) - b;

  const double error = std::fabs( value - exact_value);
  const double tol = 1.0e-15;

  int return_flag = 0;

  if( error > tol )
    {
      std::cout << "Error: GaussianXYProfile tolerance exceeded." << std::endl
                << "exact value = " << exact_value << std::endl
                << "value = " << value << std::endl
                << "error = " << error << std::endl
                << "tolerance = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}
