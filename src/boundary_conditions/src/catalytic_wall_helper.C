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
#include "grins/catalytic_wall_helper.h"

// GRINS
#include "grins/math_constants.h"

namespace GRINS
{
  CatalyticWallHelper::CatalyticWallHelper( const libMesh::Real R_s, const libMesh::Real M_s,
					    const libMesh::Real gamma_s )
    : _gamma_s(gamma_s),
      _C( std::sqrt( R_s/(GRINS::Constants::two_pi*M_s) ) )
  {
    return;
  }

  CatalyticWallHelper::~CatalyticWallHelper()
  {
    return;
  }

} // end namespace GRINS
