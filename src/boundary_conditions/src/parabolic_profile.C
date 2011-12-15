//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

GRINS::ParabolicProfile::ParabolicProfile( const GRINS::VariableIndex u_var_in )
  : DirichletFuncObj(),
    _u_var( u_var_in ),
    _a(0.0), _b(0.0), _c(-4.0), _d(0.0), _e(4.0), _f(0.0)
{
  return;
}

GRINS::ParabolicProfile::ParabolicProfile( const GRINS::VariableIndex u_var_in,
					   const double a, const double b, const double c,
					   const double d, const double e, const double f )
  : DirichletFuncObj(),
    _u_var( u_var_in ),
    _a(a), _b(b), _c(c), _d(d), _e(e), _f(f)
{
  return;
}

GRINS::ParabolicProfile::~ParabolicProfile()
{
  return;
}

libMesh::Number GRINS::ParabolicProfile::value( const libMesh::FEMContext& c, const unsigned int qp )
{
  const std::vector<libMesh::Point>& qpoint = c.side_fe_var[_u_var]->get_xyz();

  const double x = qpoint[qp](0);
  const double y = qpoint[qp](1);
  
  return this->eval( _a, _b, _c, _d, _e, _f, x, y );
}
