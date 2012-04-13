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

#include "concentric_cylinder_profile.h"

GRINS::ConcentricCylinderProfile::ConcentricCylinderProfile( )
  : FunctionBase<Number>(),
    _u0(2.0),
    _r0(1.0),
    _r1(2.0)
{
  return;
}

GRINS::ConcentricCylinderProfile::ConcentricCylinderProfile( const double u0, 
							     const double r0, 
							     const double r1 )
  : FunctionBase<Number>(),
    _u0(u0),
    _r0(r0),
    _r1(r1)
{
  return;
}

GRINS::ConcentricCylinderProfile::~ConcentricCylinderProfile()
{
  return;
}

libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > GRINS::ConcentricCylinderProfile::clone() const
{
  return libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> >( new ConcentricCylinderProfile( _u0, _r0, _r1 ) );
}

libMesh::Number GRINS::ConcentricCylinderProfile::operator()( const Point &p, 
							      const Real )
{
  const double r = p(0);
  
  return this->eval( _u0, _r0, _r1, r );
}

void GRINS::ConcentricCylinderProfile::operator()( const Point &p, 
						   const Real time, 
						   libMesh::DenseVector<libMesh::Number> &output )
{
  output(0) = (*this)(p, time);

  return;
}
