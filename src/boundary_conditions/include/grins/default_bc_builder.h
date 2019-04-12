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

#ifndef GRINS_DEFAULT_BC_BUILDER_H
#define GRINS_DEFAULT_BC_BUILDER_H

// GRINS
#include "grins/bc_builder.h"
#include "grins/builder_helper.h"

// libMesh foward declarations
namespace libMesh
{
  class System;
}

namespace GRINS
{
  //! Manages runtime construction of Dirichlet boundary conditions
  /*! This will parse the input for the request Dirichlet boundary
    conditions and manage their construction. Actual construction of
    the DirichletBoundary objects is delegated to factory
    classes. This builder classes merely manages tasks around the
    factories as needed.  To add new Dirichlet boundary conditions,
    the user should instantiate an appropriate factory sub class. */
  class DefaultBCBuilder : public BCBuilder,
                           public BuilderHelper
  {
  public:

    DefaultBCBuilder()
      : BCBuilder(),
        BuilderHelper()
    {};

    ~DefaultBCBuilder(){};

  protected:

    virtual void build_bcs( const GetPot& input, MultiphysicsSystem& system,
                            std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

    //! Helper function to build boundary conditions specified by a single type
    /*! Examples include axisymmetric and periodic. */
    void build_type_based_bcs( const GetPot& input,
                               MultiphysicsSystem& system,
                               const std::set<BoundaryID>& bc_ids,
                               libMesh::DofMap& dof_map,
                               const std::string& type_input_section,
                               std::set<std::string>& var_sections,
                               std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

    void build_axisymmetric_bcs( const GetPot& input,
                                 MultiphysicsSystem& system,
                                 const std::set<BoundaryID>& bc_ids,
                                 libMesh::DofMap& dof_map,
                                 const std::string& bc_type,
                                 std::set<std::string>& var_sections,
                                 std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

    //! Helper function to build boundary conditions using Variable sections
    /*! This is the "standard" part. We parse for each Variable section that
      should have boundary conditions and then parse the boundary condition
      type. */
    void build_bcs_by_var_section(const GetPot& input,
                                  MultiphysicsSystem& system,
                                  const std::string& bc_name,
                                  const std::set<BoundaryID>& bc_ids,
                                  libMesh::DofMap& dof_map,
                                  std::set<std::string>& var_sections,
                                  const std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >& bc_id_to_subdomain_id_map,
                                  std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs);

    void parse_and_build_bc_id_map( const GetPot& input,
                                    std::map<std::string,std::set<BoundaryID> >& bc_id_map );

    void verify_bc_ids_with_mesh( const MultiphysicsSystem& system,
                                  const std::map<std::string,std::set<BoundaryID> >& bc_id_map ) const;

    void build_periodic_bc( const GetPot& input,
                            libMesh::System& system,
                            const std::set<BoundaryID>& bc_ids,
                            const std::string& section );

    void parse_periodic_master_slave_ids( const GetPot& input,
                                          const std::string& section,
                                          libMesh::boundary_id_type& master_id,
                                          libMesh::boundary_id_type& slave_id ) const;

    libMesh::RealVectorValue parse_periodic_offset(const GetPot& input,
                                                   const std::string& section) const;

    //! Build up bc_id to subdomain_id map
    /*! we also check and make sure that there's only one subdomain id
      per boundary id. If not, we throw an error. */
    void build_bc_to_subdomain_map_check_with_mesh
    ( const MultiphysicsSystem& system,
      std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >& bc_id_to_subdomain_id_map ) const;

    //! Check if the Variable var is active on the given subdomain_id
    bool is_var_active( const FEVariablesBase& var,
                        const std::vector<libMesh::subdomain_id_type>& subdomain_ids ) const;

  };
} // end namespace GRINS

#endif // GRINS_DEFAULT_BC_BUILDER_H
