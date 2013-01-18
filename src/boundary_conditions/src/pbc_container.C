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

#include "grins/pbc_container.h"

namespace GRINS
{

  PBCContainer::PBCContainer()
    : _master_id( -1 ),
      _slave_id( -1 ),
      _offset_vector( 0.0, 0.0, 0.0 )
  {
    return;
  }

  PBCContainer::~PBCContainer()
  {
    return;
  }

  void PBCContainer::set_master_bcid( const BoundaryID bc_id )
  {
    _master_id = bc_id;
    return;
  }

  void PBCContainer::set_slave_bcid( const BoundaryID bc_id )
  {
    _slave_id = bc_id;
    return;
  }

  void PBCContainer::set_offset_vector( const libMesh::RealVectorValue& offset_vector )
  {
    _offset_vector = offset_vector;
    return;
  }

  BoundaryID PBCContainer::get_master_bcid() const
  {
    return _master_id;
  }

  BoundaryID PBCContainer::get_slave_bcid() const
  {
    return _slave_id;
  }

  const libMesh::RealVectorValue& PBCContainer::get_offset_vector() const
  {
    return _offset_vector;
  }

} // namespace GRINS
