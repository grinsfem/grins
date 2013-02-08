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
#ifndef GRINS_BC_HANDLING_BASE_H
#define GRINS_BC_HANDLING_BASE_H

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/point.h"

//GRINS
#include "grins/variable_name_defaults.h"
#include "grins/var_typedefs.h"
#include "grins/bc_types.h"
#include "grins/dbc_container.h"
#include "grins/pbc_container.h"
#include "grins/nbc_container.h"
#include "grins/boundary_conditions.h"
#include "grins/cached_values.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMContext;
  class FEMSystem;
  class DofMap;
}

namespace GRINS
{
  //! Base class for reading and handling boundary conditions for physics classes
  class BCHandlingBase
  {
  public:
    
    BCHandlingBase(const std::string& physics_name);
    
    virtual ~BCHandlingBase();

    void attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs );

    void attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc );

    virtual void read_bc_data( const GetPot& input, const std::string& id_str,
			       const std::string& bc_str );

    //! Override this method to initialize any system-dependent data.
    /*! Override this method to, for example, cache a System variable
        number. This is called before any of the other init methods in this class.
        By default, does nothing. */
    virtual void init_bc_data( const libMesh::FEMSystem& system );

    virtual void apply_neumann_bcs( libMesh::FEMContext& context,
				    const GRINS::CachedValues& cache,
				    const bool request_jacobian,
				    const GRINS::BoundaryID bc_id ) const;

    virtual void user_apply_neumann_bcs( libMesh::FEMContext& context,
					 const GRINS::CachedValues& cache,
					 const bool request_jacobian,
					 const GRINS::BoundaryID bc_id,
					 const GRINS::BCType bc_type ) const;

    virtual void init_dirichlet_bc_func_objs( libMesh::FEMSystem* system ) const;

    virtual void init_periodic_bcs( libMesh::FEMSystem* system ) const;

    void set_dirichlet_bc_type( GRINS::BoundaryID bc_id, int bc_type );
    void set_neumann_bc_type( GRINS::BoundaryID bc_id, int bc_type );
    void set_dirichlet_bc_value( GRINS::BoundaryID bc_id, libMesh::Real value, int component = 0 );
    void set_neumann_bc_value( GRINS::BoundaryID bc_id, const libMesh::Point& q_in );
    libMesh::Real get_dirichlet_bc_value( GRINS::BoundaryID bc_id, int component = 0 ) const;

    const libMesh::Point& get_neumann_bc_value( GRINS::BoundaryID bc_id ) const;

    std::tr1::shared_ptr< GRINS::NeumannFuncObj > get_neumann_bound_func( GRINS::BoundaryID bc_id,
									  GRINS::VariableIndex var_id ) const;

    std::tr1::shared_ptr< GRINS::NeumannFuncObj > get_neumann_bound_func( GRINS::BoundaryID bc_id,
									  GRINS::VariableIndex var_id );

    virtual void init_dirichlet_bcs( libMesh::FEMSystem* system ) const;

    // User will need to implement these functions for BC handling
    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_types( const GRINS::BoundaryID bc_id, 
				const std::string& bc_id_string, 
				const int bc_type, 
				const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system,
					  libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id,
					  GRINS::BCType bc_type ) const;
    

    GRINS::BCType get_dirichlet_bc_type( const GRINS::BoundaryID bc_id ) const;

    bool is_axisymmetric() const;

  protected:

    //! Map between boundary id and Dirichlet boundary condition type
    /*! We need to keep this around because the libMesh::DirichletBoundary
      objects can't be created until we init the variables */
    std::map< GRINS::BoundaryID, GRINS::BCType> _dirichlet_bc_map;

    //! Map between boundary id and Neumann boundary condition type
    std::map< GRINS::BoundaryID, GRINS::BCType> _neumann_bc_map;

    //! Stash prescribed Dirichlet boundary values
    std::map< GRINS::BoundaryID, libMesh::Point > _dirichlet_values;

    //! Stash prescribed boundary fluxes
    std::map< GRINS::BoundaryID, libMesh::Point > _q_values;

    
    std::map< GRINS::BoundaryID, GRINS::NBCContainer > _neumann_bound_funcs;

    std::vector< GRINS::DBCContainer > _dirichlet_bound_funcs;

    //! Object that stashes generic boundary condition types
    /** \todo Move this so that only one object is needed. 
	      Perhaps make static? */
    GRINS::BoundaryConditions _bound_conds;

    std::vector< GRINS::PBCContainer > _periodic_bcs;
    unsigned int _num_periodic_bcs;

    std::string _physics_name;

    enum BC_BASE{ PERIODIC = -1};

    //! Flag to cache whether or not there is an axisymmetric boundary present
    static bool _axisymmetric;

  private:
    BCHandlingBase();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  bool BCHandlingBase::is_axisymmetric() const
  {
    return _axisymmetric;
  }

  inline
  const libMesh::Point& BCHandlingBase::get_neumann_bc_value( GRINS::BoundaryID bc_id ) const
  {
    return (_q_values.find(bc_id))->second;
  }

  inline
  std::tr1::shared_ptr< GRINS::NeumannFuncObj >
  BCHandlingBase::get_neumann_bound_func( GRINS::BoundaryID bc_id,
					  GRINS::VariableIndex var_id ) const
  {
    return ((_neumann_bound_funcs.find(bc_id))->second).get_func(var_id);
  }

  inline
  std::tr1::shared_ptr< GRINS::NeumannFuncObj >
  BCHandlingBase::get_neumann_bound_func( GRINS::BoundaryID bc_id,
					  GRINS::VariableIndex var_id )
  {
    return ((_neumann_bound_funcs.find(bc_id))->second).get_func(var_id);
  }
  
  inline
  GRINS::BCType BCHandlingBase::get_dirichlet_bc_type( const GRINS::BoundaryID bc_id ) const
  {
    std::map< GRINS::BoundaryID, GRINS::BCType>::const_iterator it = 
      _dirichlet_bc_map.find(bc_id);
    return it->second;
  }

}
#endif // GRINS_BC_HANDLING_BASE
