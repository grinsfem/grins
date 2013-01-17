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

#include "grins/gaussian_xy_profile.h"

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

  libMesh::Number GRINS::GaussianXYProfile::operator()( const Point &p, 
							const Real )
  {
    const double r = std::sqrt( p(0)*p(0) + p(1)*p(1) );
  
    return this->eval( r );
  }

  void GRINS::GaussianXYProfile::operator()( const Point &p, 
					     const Real time, 
					     libMesh::DenseVector<libMesh::Number> &output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p, time);
      }
    return;
  }

  libMesh::Number GRINS::GaussianXYProfile::operator()( unsigned int i,
							const Point &p, 
							const Real time )
  {
    return (*this)(p, time);
  }

} // namespace GRINS
