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

#ifndef GRINS_REGRESSION_HELPER_H
#define GRINS_REGRESSION_HELPER_H

// C++
#include <limits>

// libMesh
#include "libmesh/libmesh_common.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>

#endif // GRINS_HAVE_CPPUNIT


namespace GRINSTesting
{
  //! Provides helper functions for fitting data to a curve.
  //!
  //! Currently supports linear regression \f$ Y = m X + b \f$
  class RegressionHelper
  {
  protected:

    //! Slope \f$ m \f$ for linear regression
    libMesh::Real linear_slope(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y, unsigned int start_index, unsigned int end_index)
    {
      libmesh_assert(x.size() == y.size());
      libmesh_assert(start_index < end_index);

      unsigned int n = x.size();
      libmesh_assert(end_index < n);

      libMesh::Real sum_y = 0.0,
        sum_x2 = 0.0,
        sum_x = 0.0,
        sum_xy = 0.0;

      for (unsigned int i=start_index; i<=end_index; i++)
        {
          libMesh::Real xi = x[i];
          libMesh::Real yi = y[i];

          sum_y  += yi;
          sum_x2 += xi*xi;
          sum_x  += xi;
          sum_xy += xi*yi;
        }


      libMesh::Real slope = ( n*sum_xy - sum_x*sum_y )/( n*sum_x2 - (sum_x*sum_x) );
      return slope;
    }

    //! y-intercept \f$ b \f$ for linear regression
    libMesh::Real linear_intercept(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y, unsigned int start_index, unsigned int end_index)
    {
      libmesh_assert(x.size() == y.size());
      libmesh_assert(start_index < end_index);

      unsigned int n = x.size();
      libmesh_assert(end_index < n);

      libMesh::Real sum_y = 0.0,
        sum_x2 = 0.0,
        sum_x = 0.0,
        sum_xy = 0.0;

      for (unsigned int i=start_index; i<=end_index; i++)
        {
          libMesh::Real xi = x[i];
          libMesh::Real yi = y[i];

          sum_y  += yi;
          sum_x2 += xi*xi;
          sum_x  += xi;
          sum_xy += xi*yi;
        }

      libMesh::Real intercept = ( sum_y*sum_x2 - sum_x*sum_xy )/( n*sum_x2 - (sum_x*sum_x) );
      return intercept;
    }

    //! Check the convergence rate of the supplied data to within the given absolute tolerance
    void check_convergence_rate(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y, unsigned int start_index, unsigned int end_index, libMesh::Real convergence_rate, libMesh::Real tol)
    {
      libmesh_assert(tol > 0.0);
      libmesh_assert(x.size() == y.size());
      libmesh_assert(start_index < end_index);
      libmesh_assert(end_index < x.size());

      libMesh::Real calculated_rate = linear_slope(x,y,start_index,end_index);

#ifdef GRINS_HAVE_CPPUNIT

      CPPUNIT_ASSERT_DOUBLES_EQUAL(convergence_rate,calculated_rate,tol);
#else
      if (std::abs(convergence_rate - calculated_rate) > tol)
        {
          std::stringstream ss;
          ss <<"ERROR" <<std::endl
             <<"calculated convergence rate: " <<calculated_rate <<std::endl
             <<"desired convergence rate: " <<convergence_rate <<std::endl
             <<"tolerance: " <<tol <<std::endl;

          libmesh_error_msg(ss.str());
        }
#endif // GRINS_HAVE_CPPUNIT
    }

  };
}

#endif // GRINS_REGRESSION_HELPER_H
