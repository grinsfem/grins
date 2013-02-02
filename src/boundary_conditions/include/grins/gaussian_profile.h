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
#ifndef GAUSSIAN_PROFILE_H
#define GAUSSIAN_PROFILE_H

// GRINS
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/function_base.h"

namespace GRINS
{
  
  //! Gaussian profile
  /*! Generic Gaussian profile. Mainly used for defining inflow
      boundary conditions. Function takes the form:
      \f$ f= a \exp \left\{ -\frac{(r-\mu)^2}{2*\sigma^2} \right\} - b\f$ 
      where \f$ b \f$ is a shift so that \f$ f \f$ is continuous.
      In particular, if \f$ r > b, f = 0 \f$.  \f$ r \f$ is computed
      in the derived classes, depending on the plane the function is defined in.*/
  
  class GaussianProfile : public libMesh::FunctionBase<libMesh::Number>
  {
  public:

    GaussianProfile( const double a, const double mu, const double sigma,
		     const double b );

    virtual ~GaussianProfile( );
    
  protected:
    
    inline double eval( const double r )
    {
      const double r_m_mu = r - _mu;
      return _a*std::exp( -r_m_mu*r_m_mu/(2.0*_variance) ) - _b;
    };

    //! Coefficients defining parabola
    double _a, _mu, _variance, _b;

  private:
    
    GaussianProfile();

  }; // class GaussianProfile

} // namespace GRINS

#endif // GAUSSIAN_PROFILE_H
