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

#ifndef GRINS_OLD_STYLE_BC_BUILDER_H
#define GRINS_OLD_STYLE_BC_BUILDER_H

// GRINS
#include "grins/bc_builder.h"

namespace GRINS
{
  //! Manages runtime construction of Dirichlet boundary conditions
  /*! This will parse the input for the request Dirichlet boundary
    conditions and manage their construction. Actual construction of
    the DirichletBoundary objects is delegated to factory
    classes. This builder classes merely manages tasks around the
    factories as needed.  To add new Dirichlet boundary conditions,
    the user should instantiate an appropriate factory sub class. */
  class OldStyleBCBuilder : public BCBuilder
  {
  public:

    OldStyleBCBuilder()
      : BCBuilder()
    {};

    ~OldStyleBCBuilder(){};

  protected:

    virtual void build_bcs( const GetPot& input, MultiphysicsSystem& system,
                            std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

  private:

    //! Determine the FEVariable type from the raw_physics_name
    /*! This will set var_section and return the FEVariablesBase pointer. We do
      it this way, instead of a void function, to avoid a [-Wunused-but-set-parameter]
      warning on setting the FEVariablesBase pointer. */
    const FEVariablesBase* determine_variable_group( const std::string& raw_physics_name,
                                                     const std::string& bc_type_str,
                                                     std::string& var_section );

    void construct_bcs_old_style( const GetPot& input,
                                  MultiphysicsSystem& system,
                                  const std::string& raw_physics_name,
                                  const std::string& section_name,
                                  const std::string& bc_id_str,
                                  const std::string& bc_type_str,
                                  const std::string& bc_value_str,
                                  const std::string& bc_var_str,
                                  libMesh::DofMap& dof_map,
                                  std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs );

    void build_basic_physics( std::set<std::string>& physics_names );

    void build_vel_and_temp_physics( std::set<std::string>& physics_names );

    void build_reacting_physics( std::set<std::string>& physics_names );

    template<typename FunctionType>
    void set_dirichlet_bc_factory_old_style_quantities( const std::string& bc_value_str,
                                                        unsigned int value_idx,
                                                        const std::vector<std::string>& var_names );

    template<typename FunctionType>
    void set_neumann_bc_factory_old_style_quantities( const std::string& bc_value_str,
                                                      unsigned int value_idx,
                                                      const std::vector<std::string>& var_names );

    void build_periodic_bc( const GetPot& input,
                            const std::string& section_name,
                            BoundaryID bc_id,
                            libMesh::DofMap& dof_map );

  };
} // end namespace GRINS

#endif // GRINS_OLD_STYLE_BC_BUILDER_H
