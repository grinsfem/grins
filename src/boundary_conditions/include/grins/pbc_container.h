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

#ifndef PBC_CONTAINER_H
#define PBC_CONTAINER_H

// libMesh
#include "libmesh/vector_value.h"

// GRINS
#include "grins/var_typedefs.h"

namespace GRINS
{
  //! Simple helper class to setup periodic boundary conditions
  /*! This class is to temporarily stash data necessary for setting
      up libMesh::PeriodicBoundary objects. Actually instantiation
      of libMesh::DirichletBoundary object is handled internally by
      GRINS::BCHandlingBase::init_periodic_bcs. For each periodic bc pair
      there is a unique PBCContainer object. */
  class PBCContainer
  {
  public:
    
    PBCContainer();
    ~PBCContainer();

    //! Add variables that are constrained by the Dirichlet bc
    void set_master_bcid( const GRINS::BoundaryID bc_id );
    void set_slave_bcid( const GRINS::BoundaryID bc_id );

    void set_offset_vector( const libMesh::RealVectorValue& offset_vector );

    GRINS::BoundaryID get_master_bcid() const;
    GRINS::BoundaryID get_slave_bcid() const;
    const libMesh::RealVectorValue& get_offset_vector() const;

  private:
    
    GRINS::BoundaryID _master_id, _slave_id;
    libMesh::RealVectorValue _offset_vector;

  };
}
#endif //PBC_CONTAINER_H
