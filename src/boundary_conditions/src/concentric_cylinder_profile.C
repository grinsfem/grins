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
#include "grins/concentric_cylinder_profile.h"

// libMesh
#include "libmesh/point.h"

namespace GRINS
{

  ConcentricCylinderProfile::ConcentricCylinderProfile( )
    : libMesh::FunctionBase<libMesh::Number>(),
      _u0(2.0),
      _r0(1.0),
      _r1(2.0)
  {
    return;
  }

  ConcentricCylinderProfile::ConcentricCylinderProfile( const double u0, 
							const double r0, 
							const double r1 )
    : libMesh::FunctionBase<libMesh::Number>(),
      _u0(u0),
      _r0(r0),
      _r1(r1)
  {
    return;
  }

  ConcentricCylinderProfile::~ConcentricCylinderProfile()
  {
    return;
  }

  libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > ConcentricCylinderProfile::clone() const
  {
    return libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> >( new ConcentricCylinderProfile( _u0, _r0, _r1 ) );
  }

  libMesh::Number ConcentricCylinderProfile::operator()( const libMesh::Point& p, 
							 const libMesh::Real )
  {
    const double r = p(0);
  
    return this->eval( _u0, _r0, _r1, r );
  }

  void ConcentricCylinderProfile::operator()( const libMesh::Point& p, 
					      const libMesh::Real time, 
					      libMesh::DenseVector<libMesh::Number> &output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p, time);
      }
    return;
  }

  libMesh::Number ConcentricCylinderProfile::component( unsigned int /*i*/,
							const libMesh::Point& p, 
							const libMesh::Real time )
  {
    return (*this)(p, time);
  }

} // namespace GRINS
