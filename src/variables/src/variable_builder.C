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

// This class
#include "grins/variable_builder.h"

// GRINS
#include "grins/default_variable_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/variable_factory.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

namespace GRINS
{
  void VariableBuilder::build_variables( const GetPot& input,
                                         MultiphysicsSystem& system )
  {
    DefaultVariableBuilder var_builder;

    // Defer the construction to the builder subclass implementation
    var_builder.build_variables_impl(input,system);
  }

  void VariableBuilder::add_variable_to_warehouse( std::shared_ptr<FEVariablesBase>& fe_var,
                                                   const std::string& var_name )
  {
    GRINSPrivate::VariableWarehouse::register_variable(var_name,fe_var);
  }

  void VariableBuilder::add_vars_to_system( MultiphysicsSystem& system,
                                            const std::vector<std::string>& var_names,
                                            const std::string& fe_family,
                                            const std::string& order,
                                            std::vector<VariableIndex>& var_indices,
                                            const std::set<libMesh::subdomain_id_type>& subdomain_ids )
  {
    const unsigned int n_vars = var_names.size();

    // Setup var_indices
    libmesh_assert( var_indices.empty() );
    var_indices.resize(n_vars);

    if( subdomain_ids.empty() )
      for( unsigned int v = 0; v < n_vars; v++ )
        var_indices[v] = system.add_variable( var_names[v],
                                              libMesh::Utility::string_to_enum<GRINSEnums::Order>(order),
                                              libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(fe_family) );

    else
      for( unsigned int v = 0; v < n_vars; v++ )
        var_indices[v] = system.add_variable( var_names[v],
                                              libMesh::Utility::string_to_enum<GRINSEnums::Order>(order),
                                              libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(fe_family),
                                              &subdomain_ids );
  }

  std::shared_ptr<FEVariablesBase> VariableBuilder::build_fe_var( const std::string& var_type,
                                                            const std::vector<std::string>& var_names,
                                                            const std::vector<VariableIndex>&  var_indices,
                                                            const std::set<libMesh::subdomain_id_type>& subdomain_ids )
  {
    // Setup VariableFactory
    VariableFactoryAbstract::set_var_names(var_names);
    VariableFactoryAbstract::set_var_indices(var_indices);
    VariableFactoryAbstract::set_subdomain_ids(subdomain_ids);

    std::unique_ptr<FEVariablesBase> var = VariableFactoryAbstract::build(var_type);

    // Need to return a std::shared_ptr, so release from the UniquePtr we got back
    return std::shared_ptr<FEVariablesBase>( var.release() );
  }
} // end namespace GRINS
