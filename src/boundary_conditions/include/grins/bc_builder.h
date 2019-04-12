//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_BC_BUILDER_H
#define GRINS_BC_BUILDER_H

// C++
#include <set>
#include <string>
#include <vector>

// GRINS
#include "grins/var_typedefs.h"
#include "grins/neumann_bc_container.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/hp_coarsentest.h" // RealVectorValue

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class DofMap;
}

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  //! Manages runtime construction of Dirichlet boundary conditions
  /*! This will parse the input for the request Dirichlet boundary
    conditions and manage their construction. Actual construction of
    the DirichletBoundary objects is delegated to factory
    classes. This builder classes merely manages tasks around the
    factories as needed.  To add new Dirichlet boundary conditions,
    the user should instantiate an appropriate factory sub class. */
  class BCBuilder
  {
  public:

    BCBuilder(){};

    virtual ~BCBuilder(){};

    static void build_boundary_conditions( const GetPot& input,
                                           MultiphysicsSystem& system,
                                           std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

  protected:

    virtual void build_bcs( const GetPot& input, MultiphysicsSystem& system,
                            std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs ) =0;

    void construct_dbc_core( const GetPot& input,
                             MultiphysicsSystem& system,
                             const std::set<BoundaryID>& bc_ids,
                             const FEVariablesBase& fe_var,
                             const std::string& section,
                             const std::string& bc_type,
                             libMesh::DofMap& dof_map );

    void construct_nbc_core( const GetPot& input,
                             MultiphysicsSystem& system,
                             const std::set<BoundaryID>& bc_ids,
                             const FEVariablesBase& fe_var,
                             const std::string& section,
                             const std::string& bc_type,
                             std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

    static bool is_new_bc_input_style( const GetPot& input );

    static std::unique_ptr<BCBuilder> build_builder( const GetPot& input );

    bool is_dirichlet_bc_type( const std::string& bc_type );

    bool is_neumann_bc_type( const std::string& bc_type );

    void add_periodic_bc_to_dofmap( libMesh::boundary_id_type master_id,
                                    libMesh::boundary_id_type slave_id,
                                    const libMesh::RealVectorValue& offset_vector,
                                    libMesh::DofMap& dof_map );



  };
} // end namespace GRINS

#endif // GRINS_BC_BUILDER_H
