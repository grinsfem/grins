//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_VARIABLE_BUILDER_H
#define GRINS_VARIABLE_BUILDER_H

// GRINS
#include <memory>
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/auto_ptr.h" // std::unique_ptr

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  //! Manages runtime construction of the FEVariableBase objects
  /*! build_variables() is intended to be called before building the
    Physics and boundary conditions so that the VariableWarhouse is
    populated and can be referenced by the Physics and boundary conditions. */
  class VariableBuilder
  {
  public:

    VariableBuilder(){}
    ~VariableBuilder(){}

    static void build_variables( const GetPot& input,
                                 MultiphysicsSystem& system );

    //! Implementation of Variable construction done in subclasses
    virtual void build_variables_impl( const GetPot& input,
                                       MultiphysicsSystem& system ) =0;

  protected:

    //! Adds/registers the fe_var to VariableWarehouse
    void add_variable_to_warehouse( std::shared_ptr<FEVariablesBase>& fe_var,
                                    const std::string& var_name );

    //! Given the names, family, and order, this adds the variables to the system and populates var_indices
    /*! The var_indices are the respective indices returned by the System from the add_variable
      call. */
    void add_vars_to_system( MultiphysicsSystem& system,
                             const std::vector<std::string>& var_names,
                             const std::string& fe_family,
                             const std::string& order,
                             std::vector<VariableIndex>& var_indices,
                             const std::set<libMesh::subdomain_id_type>& subdomain_ids);


    //! Sets appropriate data in the VariableFactoryAbstract and calls VariableFactoryAbstract::build()
    std::shared_ptr<FEVariablesBase> build_fe_var( const std::string& var_type,
                                             const std::vector<std::string>& var_names,
                                             const std::vector<VariableIndex>&  var_indices,
                                             const std::set<libMesh::subdomain_id_type>& subdomain_ids );

  };
} // end namespace GRINS

#endif // GRINS_VARIABLE_BUILDER_H
