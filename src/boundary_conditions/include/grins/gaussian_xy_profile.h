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
#ifndef GAUSSIAN_XY_PROFILE_H
#define GAUSSIAN_XY_PROFILE_H

// GRINS
#include "grins/gaussian_profile.h"

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

    virtual libMesh::Number operator()( const libMesh::Point &p, const libMesh::Real time );

    virtual void operator()( const libMesh::Point &p, 
			     const libMesh::Real time, 
			     libMesh::DenseVector<libMesh::Number> &output );

    virtual libMesh::Number component( unsigned int i, const libMesh::Point &p, 
				       const libMesh::Real time );

  private:
    
    GaussianXYProfile();

  }; // class GaussianXYProfile

} // namespace GRINS

#endif // GAUSSIAN_XY_PROFILE_H
