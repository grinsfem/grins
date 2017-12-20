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
#include "grins/default_variable_builder.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/variable_factory.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh_base.h"

namespace GRINS
{
  void DefaultVariableBuilder::build_variables_impl( const GetPot& input,
                                                     MultiphysicsSystem& system )
  {
    std::set<libMesh::subdomain_id_type> mesh_subdomain_ids;
    const libMesh::MeshBase& mesh = system.get_mesh();

    // Note this does a full loop of the mesh, these are not cached
    // So we compute them once
    mesh.subdomain_ids(mesh_subdomain_ids);

    std::vector<std::string> var_sections;
    this->parse_var_sections_vector( input, var_sections );

    for( std::vector<std::string>::const_iterator var_sect = var_sections.begin();
         var_sect != var_sections.end(); ++var_sect )
      {
        // Convenience
        std::string var_section( *var_sect );

        // We construct based on the Variable type so we can name the Variables whatever we want
        // If var_type is invalid, that will be detected in the build_fe_var call, since
        // we feed the var_type to VariableFactoryAbstract::build.
        std::string var_type = this->parse_var_type( input, var_section );

        // Parse names
        std::vector<std::string> var_names;
        this->parse_var_names( input, var_type, var_section, var_names );

        // Parse FE family
        std::string fe_family = this->parse_fe_family( input, var_section, var_type );

        // Parse FE order
        std::string order = this->parse_fe_order( input, var_section, var_type );

        // Parse the subdomain_ids
        std::set<libMesh::subdomain_id_type> subdomain_ids;
        this->parse_subdomain_ids( mesh_subdomain_ids, input, var_section, subdomain_ids );

        // Add variables to system
        std::vector<VariableIndex> var_indices;
        this->add_vars_to_system( system, var_names, fe_family, order, var_indices, subdomain_ids );


        // Build FEVariablesBase object
        std::shared_ptr<FEVariablesBase> fe_var = this->build_fe_var( var_type, var_names, var_indices, subdomain_ids );

        // Add to VariableWarehouse
        this->add_variable_to_warehouse( fe_var, var_section );
      }
  }

  void DefaultVariableBuilder::parse_var_names( const GetPot& input,
                                                const std::string& var_type,
                                                const std::string& var_section,
                                                std::vector<std::string>& var_names ) const
  {
    // Just in case
    var_names.clear();

    VariableFactoryAbstract::set_getpot(input);
    VariableFactoryAbstract::set_var_section(VariablesParsing::variables_section()+"/"+var_section);

    var_names = VariableFactoryAbstract::build_var_names(var_type);
  }

  std::string DefaultVariableBuilder::parse_fe_family( const GetPot& input,
                                                       const std::string& var_section,
                                                       const std::string& var_type ) const
  {
    VariableFactoryAbstract::set_getpot(input);
    VariableFactoryAbstract::set_var_section(VariablesParsing::variables_section()+"/"+var_section);

    return VariableFactoryAbstract::parse_fe_family(var_type);
  }

  std::string DefaultVariableBuilder::parse_fe_order( const GetPot& input,
                                                      const std::string& var_section,
                                                      const std::string& var_type ) const
  {
    VariableFactoryAbstract::set_getpot(input);
    VariableFactoryAbstract::set_var_section(VariablesParsing::variables_section()+"/"+var_section);

    return VariableFactoryAbstract::parse_fe_order(var_type);
  }

  void DefaultVariableBuilder::parse_subdomain_ids( const std::set<libMesh::subdomain_id_type>& mesh_subdomain_ids,
                                                    const GetPot& input,
                                                    const std::string& var_section,
                                                    std::set<libMesh::subdomain_id_type>& subdomain_ids )
  {
    std::string input_sec(VariablesParsing::variables_section()+"/"+var_section+"/enabled_subdomains");
    if( input.have_variable(input_sec) )
      {
        unsigned int invalid_id = std::numeric_limits<unsigned int>::max();

        unsigned int n_ids = input.vector_variable_size(input_sec);
        for( unsigned int i = 0; i < n_ids; i++ )
          {
            unsigned int id = input(input_sec,invalid_id,i);
            subdomain_ids.insert(id);
          }

        for( std::set<libMesh::subdomain_id_type>::const_iterator it = subdomain_ids.begin(); it != subdomain_ids.end(); ++it )
          if( mesh_subdomain_ids.find(*it) == mesh_subdomain_ids.end() )
            libmesh_error_msg("ERROR: Could not find subdomain " << *it << " in the mesh!");
      }
  }

} // end namespace GRINS
