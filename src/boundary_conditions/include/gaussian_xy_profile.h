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
#ifndef GAUSSIAN_XY_PROFILE_H
#define GAUSSIAN_XY_PROFILE_H

// GRINS
#include "gaussian_profile.h"

namespace GRINS
{
  
  //! Gaussian profile
  /*! \f$ r \f$ is computed in the x-y plane: \f$ r = \sqrt{x^2 + y^2} \f$ */
  class GaussianXYProfile : public GaussianProfile
  {
  public:

    GaussianXYProfile( const double a, const double mu, const double sigma,
		       const double b );

    virtual ~GaussianXYProfile( );

    virtual libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > clone() const;

    virtual libMesh::Number operator()( const Point &p, const Real time );

    virtual void operator()( const Point &p, 
			     const Real time, 
			     libMesh::DenseVector<Number> &output );

    virtual libMesh::Number operator()( unsigned int i, const Point &p, 
					const Real time );

  private:
    
    GaussianXYProfile();

  }; // class GaussianXYProfile

} // namespace GRINS

#endif // GAUSSIAN_XY_PROFILE_H
