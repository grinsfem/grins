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

// This class
#include "grins/gaussian_xy_profile.h"

// libMesh
#include "libmesh/point.h"

namespace GRINS
{

  GRINS::GaussianXYProfile::GaussianXYProfile( const double a, const double mu, const double sigma,
					       const double b )
    : GaussianProfile(a,mu,sigma,b)
  {
    _initialized = true;
    return;
  }

  GRINS::GaussianXYProfile::~GaussianXYProfile()
  {
    return;
  }

  libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > GRINS::GaussianXYProfile::clone() const
  {
    return libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> >( new GaussianXYProfile( _a, _mu, std::sqrt(_variance), _b ) );
  }

  libMesh::Number GRINS::GaussianXYProfile::operator()( const libMesh::Point &p, 
							const libMesh::Real )
  {
    const double r = std::sqrt( p(0)*p(0) + p(1)*p(1) );
  
    return this->eval( r );
  }

  void GRINS::GaussianXYProfile::operator()( const libMesh::Point &p, 
					     const libMesh::Real time, 
					     libMesh::DenseVector<libMesh::Number> &output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p, time);
      }
    return;
  }

  libMesh::Number GRINS::GaussianXYProfile::component( unsigned int /*i*/,
						       const libMesh::Point &p, 
						       const libMesh::Real time )
  {
    return (*this)(p, time);
  }

} // namespace GRINS
