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

// This class
#include "grins/old_style_bc_builder.h"

// GRINS
#include "grins/common.h"
#include "grins/var_typedefs.h"
#include "grins/physics_naming.h"
#include "grins/multiphysics_sys.h"
#include "grins/dirichlet_bc_factory_abstract.h"
#include "grins/neumann_bc_factory_abstract.h"
#include "grins/parsed_function_neumann_old_style_bc_factory.h"
#include "grins/dirichlet_bc_factory_function_old_style_base.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/multicomponent_variable.h"

// libMesh
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  void OldStyleBCBuilder::build_bcs( const GetPot& input, MultiphysicsSystem& system,
                                     std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs )
  {
    // Warn about deprecation of this horrid style
    {
      std::string warning = "WARNING: Specifying boundary conditions in the\n";
      warning += "         Physics sections is DEPRECATED! Please\n";
      warning += "         update your input file to new the newer\n";
      warning += "         style. See the examples for an illustration.\n";

      grins_warning(warning);
    }

    libMesh::DofMap& dof_map = system.get_dof_map();

    const PhysicsList& physics_list = system.get_physics_list();

    std::set<std::string> basic_physics;
    this->build_basic_physics(basic_physics);

    std::set<std::string> vel_and_temp_physics;
    this->build_vel_and_temp_physics(vel_and_temp_physics);

    std::set<std::string> reacting_physics;
    this->build_reacting_physics(reacting_physics);

    for( PhysicsListIter physics_iter = physics_list.begin();
         physics_iter != physics_list.end();
         physics_iter++ )
      {
        std::string physics_name = physics_iter->first;
        std::string raw_physics_name = PhysicsNaming::extract_physics(physics_name);

        std::string section_name = "Physics/"+physics_name;

        if( basic_physics.find( raw_physics_name ) != basic_physics.end() )
          {
            this->construct_bcs_old_style(input,
                                          system,
                                          raw_physics_name,
                                          section_name,
                                          "bc_ids",
                                          "bc_types",
                                          "bc_values",
                                          "bc_variables",
                                          dof_map,
                                          neumann_bcs);
          }

        if( (vel_and_temp_physics.find( raw_physics_name ) != vel_and_temp_physics.end()) ||
            (reacting_physics.find( raw_physics_name ) != reacting_physics.end()) )
          {
            this->construct_bcs_old_style(input,
                                          system,
                                          raw_physics_name,
                                          section_name,
                                          "vel_bc_ids",
                                          "vel_bc_types",
                                          "vel_bc_values",
                                          "vel_bc_variables",
                                          dof_map,
                                          neumann_bcs);

            this->construct_bcs_old_style(input,
                                          system,
                                          raw_physics_name,
                                          section_name,
                                          "temp_bc_ids",
                                          "temp_bc_types",
                                          "temp_bc_values",
                                          "temp_bc_variables",
                                          dof_map,
                                          neumann_bcs);
          }

        if( reacting_physics.find( raw_physics_name ) != reacting_physics.end() )
          {
            this->construct_bcs_old_style(input,
                                          system,
                                          raw_physics_name,
                                          section_name,
                                          "species_bc_ids",
                                          "species_bc_types",
                                          "species_bc_values",
                                          "species_bc_variables",
                                          dof_map,
                                          neumann_bcs);
          }
      }
  }

  const FEVariablesBase* OldStyleBCBuilder::determine_variable_group( const std::string& raw_physics_name,
                                                                      const std::string& bc_type_str,
                                                                      std::string& var_section )
  {
    if( bc_type_str == std::string("bc_types") )
      {
        if( raw_physics_name == PhysicsNaming::incompressible_navier_stokes() ||
            raw_physics_name == PhysicsNaming::stokes() )
          var_section = VariablesParsing::velocity_section();
        else if( raw_physics_name == PhysicsNaming::heat_conduction() ||
                 raw_physics_name == PhysicsNaming::heat_transfer() ||
                 raw_physics_name == PhysicsNaming::axisymmetric_heat_transfer() )
          var_section = VariablesParsing::temperature_section();
        else if( raw_physics_name == PhysicsNaming::spalart_allmaras() )
          var_section = VariablesParsing::turbulence_section();
        else if( raw_physics_name == PhysicsNaming::elastic_membrane() ||
                 raw_physics_name == PhysicsNaming::elastic_cable() )
          var_section = VariablesParsing::displacement_section();
        else if( raw_physics_name == PhysicsNaming::convection_diffusion() )
          var_section = VariablesParsing::single_var_section()+":"+raw_physics_name;
        else
          libmesh_error();
      }
    else if( bc_type_str == std::string("vel_bc_types") )
      var_section = VariablesParsing::velocity_section();

    else if( bc_type_str == std::string("temp_bc_types") )
      var_section = VariablesParsing::temperature_section();

    else if( bc_type_str == std::string("species_bc_types") )
      var_section = VariablesParsing::species_mass_fractions_section();

    else
      libmesh_error();

    return &GRINSPrivate::VariableWarehouse::get_variable(var_section);
  }

  void OldStyleBCBuilder::construct_bcs_old_style( const GetPot& input,
                                                   MultiphysicsSystem& system,
                                                   const std::string& raw_physics_name,
                                                   const std::string& section_name,
                                                   const std::string& bc_id_str,
                                                   const std::string& bc_type_str,
                                                   const std::string& bc_value_str,
                                                   const std::string& bc_var_str,
                                                   libMesh::DofMap& dof_map,
                                                   std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs )
  {
    unsigned int num_ids = input.vector_variable_size(section_name+"/"+bc_id_str);
    unsigned int num_types = input.vector_variable_size(section_name+"/"+bc_type_str);

    if( num_ids != num_types )
      libmesh_error_msg("Error: Must specify equal number of boundary ids and boundary conditions");

    for( unsigned int i = 0; i < num_ids; i++ )
      {
        // Parse the bc type, add "_old_style" at the end to distinguish for deprecated construction
        std::string bc_type = input(section_name+"/"+bc_type_str, std::string("DIE!"), i );
        bc_type += "_old_style";

        BoundaryID bc_id = input(section_name+"/"+bc_id_str, -1, i );

        // If this is a periodic boundary condition, we can immediately
        // apply and move to the next one
        if( bc_type == std::string("periodic_old_style") )
          {
            this->build_periodic_bc( input, section_name, bc_id, dof_map );
            continue;
          }

        // We use the set for compatibility with the BCFactories
        std::set<BoundaryID> bc_ids;
        bc_ids.insert(bc_id);

        std::string variable_group_name;

        const FEVariablesBase* fe_var_ptr = this->determine_variable_group( raw_physics_name,
                                                                            bc_type_str,
                                                                            variable_group_name );

        libmesh_assert(fe_var_ptr);

        // We need the var_names for the old style parsing
        std::vector<std::string> var_names;
        var_names = fe_var_ptr->active_var_names();

        // Axisymmetric is special. It depends on the variable type.
        // So, we prepend the type with the variable name in that
        // case.
        if( bc_type == "axisymmetric_old_style" )
          bc_type = variable_group_name+"_"+bc_type;


        // For these types of boundary conditions, we need to treat one
        // variable at a time, so extract the relevant one.
        if( bc_type == std::string("parsed_dirichlet_old_style") ||
            bc_type == std::string("constant_dirichlet_old_style") ||
            bc_type == std::string("parsed_fem_dirichlet_old_style") ||
            bc_type == std::string("parsed_neumann_old_style") )
          {
            var_names.clear();
            var_names.resize(1, input(section_name+"/"+bc_var_str, std::string("DIE!"), i ) );
          }

        if( this->is_dirichlet_bc_type(bc_type) )
          {
            // Tell the old style DirichletBCFactory where to parse the value of the BC
            this->set_dirichlet_bc_factory_old_style_quantities<libMesh::FunctionBase<libMesh::Number> >
              ( bc_value_str, i, var_names );
            this->set_dirichlet_bc_factory_old_style_quantities<libMesh::FEMFunctionBase<libMesh::Number> >
              ( bc_value_str, i, var_names );

            this->construct_dbc_core( input, system, bc_ids, *fe_var_ptr,
                                      section_name, bc_type, dof_map );
          }
        else if( this->is_neumann_bc_type(bc_type) )
          {
            // Tell the old style NeumannBCFactory where to parse the value of the BC
            this->set_neumann_bc_factory_old_style_quantities<libMesh::FunctionBase<libMesh::Number> >
              ( bc_value_str, i, var_names );
            this->set_neumann_bc_factory_old_style_quantities<libMesh::FEMFunctionBase<libMesh::Number> >
              ( bc_value_str, i, var_names );

            this->construct_nbc_core( input, system, bc_ids, *fe_var_ptr,
                                      section_name, bc_type, neumann_bcs );
          }
        else
          libmesh_error_msg("ERROR: Invalid bc_type "+bc_type+"!");

      }
  }

  void OldStyleBCBuilder::build_basic_physics( std::set<std::string>& physics_names )
  {
    physics_names.insert(PhysicsNaming::incompressible_navier_stokes());
    physics_names.insert(PhysicsNaming::stokes());
    physics_names.insert(PhysicsNaming::elastic_membrane());
    physics_names.insert(PhysicsNaming::elastic_cable());
    physics_names.insert(PhysicsNaming::convection_diffusion());
    physics_names.insert(PhysicsNaming::spalart_allmaras());
    physics_names.insert(PhysicsNaming::axisymmetric_heat_transfer());
    physics_names.insert(PhysicsNaming::heat_conduction());
    physics_names.insert(PhysicsNaming::heat_transfer());
  }

  void OldStyleBCBuilder::build_vel_and_temp_physics( std::set<std::string>& physics_names )
  {
    physics_names.insert(PhysicsNaming::low_mach_navier_stokes());
  }

  void OldStyleBCBuilder::build_reacting_physics( std::set<std::string>& physics_names )
  {
    physics_names.insert(PhysicsNaming::reacting_low_mach_navier_stokes());
  }

  void OldStyleBCBuilder::build_periodic_bc( const GetPot& input,
                                             const std::string& section,
                                             BoundaryID bc_id,
                                             libMesh::DofMap& dof_map )
  {
    std::string wall_input = section+"/periodic_wall_";
    wall_input += StringUtilities::T_to_string<BoundaryID>(bc_id);

    if( input.have_variable(wall_input) )
      {
        libMesh::boundary_id_type invalid_bid =
          std::numeric_limits<libMesh::boundary_id_type>::max();

        libMesh::boundary_id_type slave_id = invalid_bid;
        libMesh::boundary_id_type master_id = invalid_bid;

        if( input.vector_variable_size(wall_input) != 2 )
          libmesh_error_msg("ERROR: "+wall_input+" must have only 2 components!");

        master_id = bc_id;

        if( input(wall_input,invalid_bid,0) == bc_id )
          slave_id = input(wall_input,invalid_bid,1);
        else
          slave_id = input(wall_input,invalid_bid,0);

        std::string offset_input = section+"/periodic_offset_";
        offset_input += StringUtilities::T_to_string<BoundaryID>(bc_id);

        if( !input.have_variable(offset_input) )
          libmesh_error_msg("ERROR: Could not find "+offset_input+"!");

        unsigned int n_comps = input.vector_variable_size(offset_input);

        libMesh::Real invalid_real = std::numeric_limits<libMesh::Real>::max();

        libMesh::RealVectorValue offset_vector;
        for( unsigned int i = 0; i < n_comps; i++ )
          offset_vector(i) = input(offset_input,invalid_real,i);

        this->add_periodic_bc_to_dofmap( master_id, slave_id,
                                         offset_vector, dof_map );
      }
  }

  template<typename FunctionType>
  void OldStyleBCBuilder::set_dirichlet_bc_factory_old_style_quantities( const std::string& bc_value_str,
                                                                         unsigned int value_idx,
                                                                         const std::vector<std::string>& var_names )
  {
    DirichletBCFactoryFunctionOldStyleBase<FunctionType>::set_value_var_old_style( bc_value_str );
    DirichletBCFactoryFunctionOldStyleBase<FunctionType>::set_value_index_old_style( value_idx );
    DirichletBCFactoryFunctionOldStyleBase<FunctionType>::set_var_names_old_style( var_names );
  }

  template<typename FunctionType>
  void OldStyleBCBuilder::set_neumann_bc_factory_old_style_quantities( const std::string& bc_value_str,
                                                                       unsigned int value_idx,
                                                                       const std::vector<std::string>& /*var_names*/ )
  {
    ParsedFunctionNeumannOldStyleBCFactory<FunctionType>::set_value_var_old_style( bc_value_str );
    ParsedFunctionNeumannOldStyleBCFactory<FunctionType>::set_value_index_old_style( value_idx );
    //ParsedFunctionNeumannOldStyleBCFactory<FunctionType>::set_var_names_old_style( var_names );
  }

  // Instantiate
  template void OldStyleBCBuilder::set_dirichlet_bc_factory_old_style_quantities<libMesh::FunctionBase<libMesh::Number> >( const std::string&, unsigned int, const std::vector<std::string>& );
  template void OldStyleBCBuilder::set_dirichlet_bc_factory_old_style_quantities<libMesh::FEMFunctionBase<libMesh::Number> >( const std::string&, unsigned int, const std::vector<std::string>& );
  template void OldStyleBCBuilder::set_neumann_bc_factory_old_style_quantities<libMesh::FunctionBase<libMesh::Number> >( const std::string&, unsigned int, const std::vector<std::string>& );
  template void OldStyleBCBuilder::set_neumann_bc_factory_old_style_quantities<libMesh::FEMFunctionBase<libMesh::Number> >( const std::string&, unsigned int, const std::vector<std::string>& );

} // end namespace GRINS
