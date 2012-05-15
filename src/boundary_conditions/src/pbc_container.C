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

#include "pbc_container.h"

GRINS::PBCContainer::PBCContainer()
  : _master_id( -1 ),
    _slave_id( -1 ),
    _offset_vector( 0.0, 0.0, 0.0 )
{
  return;
}

GRINS::PBCContainer::~PBCContainer()
{
  return;
}

void GRINS::PBCContainer::set_master_bcid( const GRINS::BoundaryID bc_id )
{
  _master_id = bc_id;
  return;
}

void GRINS::PBCContainer::set_slave_bcid( const GRINS::BoundaryID bc_id )
{
  _slave_id = bc_id;
  return;
}

void GRINS::PBCContainer::set_offset_vector( const libMesh::RealVectorValue& offset_vector )
{
  _offset_vector = offset_vector;
  return;
}

GRINS::BoundaryID GRINS::PBCContainer::get_master_bcid() const
{
  return _master_id;
}

GRINS::BoundaryID GRINS::PBCContainer::get_slave_bcid() const
{
  return _slave_id;
}

const libMesh::RealVectorValue& GRINS::PBCContainer::get_offset_vector() const
{
  return _offset_vector;
}
