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
#ifndef BC_HANDLING_BASE_H
#define BC_HANDLING_BASE_H

//libMesh
#include "getpot.h"
#include "libmesh.h"
#include "fem_system.h"
#include "fem_context.h"
#include "dirichlet_boundaries.h"
#include "dof_map.h"

//GRINS
#include "variable_name_defaults.h"
#include "var_typedefs.h"
#include "bc_types.h"
#include "boundary_conditions.h"
#include "grins_physics_names.h"
#include "dbc_container.h"

namespace GRINS
{
  //! Base class for reading and handling boundary conditions for physics classes
  class BCHandlingBase
  {
  public:
    
    BCHandlingBase();
    
    virtual ~BCHandlingBase();

    void attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs );

    void attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc );

    virtual void read_bc_data( const GetPot& input, const std::string& id_str,
			       const std::string& bc_str );

    void init_user_dirichlet_bcs( libMesh::FEMSystem* system ) const;

    // User will need to implement these functions for BC handling
    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const GRINS::BoundaryID bc_id, 
			       const std::string& bc_id_string, 
			       const int bc_type, 
			       const GetPot& input );

    virtual void init_dirichlet_bcs( libMesh::DofMap& dof_map ) const;

    
  protected:

    //! Map between boundary id and Dirichlet boundary condition type
    /*! We need to keep this around because the libMesh::DirichletBoundary
      objects can't be created until we init the variables */
    std::map< GRINS::BoundaryID, GRINS::BCType> _dirichlet_bc_map;

    //! Map between boundary id and Neumann boundary condition type
    std::map< GRINS::BoundaryID, GRINS::BCType> _neumann_bc_map;

    //! Stash prescribed Dirichlet boundary values
    std::map< GRINS::BoundaryID, double > _dirichlet_values;

    //! Stash prescribed boundary fluxes
    std::map< GRINS::BoundaryID, libMesh::Point > _q_values;

    //! Map between boundary id and general Neumann boundary functions
    /*! The user may wish to set a different function for each variable in the physics class. 
        By design, the user cannot set more than 1 function per variable.*/
    GRINS::NBCContainer _neumann_bound_funcs;

    std::vector< GRINS::DBCContainer > _dirichlet_bound_funcs;
  };
}
#endif // BC_HANDLING_BASE
