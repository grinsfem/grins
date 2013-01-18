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

#include "constant_with_exp_layer.h"

namespace GRINS
{
  ConstantWithExponentialLayer::ConstantWithExponentialLayer( const libMesh::Real u0, const libMesh::Real factor,
							      const libMesh::Real r0, const libMesh::Real delta )
    : _u0(u0),
      _factor(factor),
      _r0(r0),
      _delta(delta)
  {
    return;
  }

  ConstantWithExponentialLayer::~ConstantWithExponentialLayer()
  {
    return;
  }
  
} // namespace GRINS
