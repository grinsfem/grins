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

#include "parabolic_profile.h"

GRINS::ParabolicProfile::ParabolicProfile( )
  : FunctionBase<Number>(),
    _a(0.0), _b(0.0), _c(-4.0), _d(0.0), _e(4.0), _f(0.0)
{
  _initialized=true;
  return;
}

GRINS::ParabolicProfile::ParabolicProfile( const double a, const double b, const double c,
					   const double d, const double e, const double f )
  : FunctionBase<Number>(),
    _a(a), _b(b), _c(c), _d(d), _e(e), _f(f)
{
  _initialized=true;
  return;
}

GRINS::ParabolicProfile::~ParabolicProfile()
{
  return;
}

libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > GRINS::ParabolicProfile::clone() const
{
  return libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> >( new ParabolicProfile( _a, _b, _c, _d, _e, _f ) );
}

libMesh::Number GRINS::ParabolicProfile::operator()( const Point &p, 
						     const Real )
{
  const double x = p(0);
  const double y = p(1);
  
  return this->eval( _a, _b, _c, _d, _e, _f, x, y );
}

void GRINS::ParabolicProfile::operator()( const Point &p, 
					  const Real time, 
					  libMesh::DenseVector<libMesh::Number> &output )
{
  for( unsigned int i = 0; i < output.size(); i++ )
    {
      output(i) = (*this)(p, time);
    }
  return;
}

libMesh::Number GRINS::ParabolicProfile::operator()( unsigned int i,
						     const Point &p, 
						     const Real time )
{
  return (*this)(p, time);
}
