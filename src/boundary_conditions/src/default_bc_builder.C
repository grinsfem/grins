//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/default_bc_builder.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/boundary_condition_names.h"
#include "grins/string_utils.h"
#include "grins/variables_parsing.h"
#include "grins/fe_variables_base.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"

namespace GRINS
{

  void DefaultBCBuilder::build_bcs( const GetPot& input, MultiphysicsSystem& system,
                                    std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs )
  {
    libMesh::DofMap& dof_map = system.get_dof_map();

    // Setup map between boundary section name and the boundary ids
    std::map<std::string,std::set<BoundaryID> > bc_id_map;
    this->parse_and_build_bc_id_map( input, bc_id_map );

    // The user must have mapped all boundary ids present in the mesh.
    // Additionally, they can't have specified an boundary id that doesn't
    // exist.
    this->verify_bc_ids_with_mesh( system, bc_id_map );

    // Parse variable sections. These correspond to the names of FEVariables
    // We will be looking for all variable sections in all boundary sections
    // that are not "constraint" variables (i.e. FEVariablesBase::_is_constraint)
    std::set<std::string> var_sections;
    this->parse_var_sections( input, var_sections );

    // Cache the map between boundary ids and subdomain ids
    // We map to a vector because it's possible to have an "interface"
    // boundary id that touches multiple elements with differing
    // subdomain ids
    std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >
             bc_id_to_subdomain_id_map;

    this->build_bc_to_subdomain_map_check_with_mesh( system,
                                                     bc_id_to_subdomain_id_map );

    for( std::map<std::string,std::set<BoundaryID> >::const_iterator bc_it = bc_id_map.begin();
         bc_it != bc_id_map.end(); ++bc_it )
      {
        const std::string& bc_name = bc_it->first;

        const std::set<BoundaryID>& bc_ids = bc_it->second;

        // Check for special types of boundary conditions that can be
        // specified as a single type that applies to all boundaries
        std::string type_input_section = BoundaryConditionNames::bc_section()+"/"+bc_name;
        std::string type_input = type_input_section+"/type";

        if( input.have_variable(type_input) )
          {
            this->build_type_based_bcs(input,system,bc_ids,dof_map,
                                       type_input_section,var_sections,
                                       neumann_bcs);
          }
        // Otherwise, we look for each Variable section and the
        // boundary conditions specified therein
        else
          {
            this->build_bcs_by_var_section(input,system,bc_name,bc_ids,dof_map,
                                           var_sections,bc_id_to_subdomain_id_map,
                                           neumann_bcs);
          }
      }
  }

  void DefaultBCBuilder::build_type_based_bcs( const GetPot& input,
                                               MultiphysicsSystem& system,
                                               const std::set<BoundaryID>& bc_ids,
                                               libMesh::DofMap& dof_map,
                                               const std::string& type_input_section,
                                               std::set<std::string>& var_sections,
                                               std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs)
  {
    std::string type_input = type_input_section+"/type";
    std::string bc_type = input( type_input, std::string("DIE!") );

    if( bc_type == BoundaryConditionNames::axisymmetric() )
      {
        // Check and make sure the Physics thinks it's axisymmetric, otherwise error
        if( !Physics::is_axisymmetric() )
          libmesh_error_msg("ERROR: Must specify Physics/is_axisymmetric = true for axisymmetric BC!");

        // Now build the boundary condition
        this->build_axisymmetric_bcs(input,system,bc_ids,dof_map,
                                     bc_type,var_sections,neumann_bcs);
      }
    else if( bc_type == BoundaryConditionNames::periodic() )
      {
        this->build_periodic_bc(input,system,bc_ids,type_input_section);
      }
    else
      {
        std::string error_msg = "ERROR: Invalid type '"+bc_type+"' for "+type_input+"!\n";
        error_msg += "       Valid values are: "+BoundaryConditionNames::axisymmetric()+"\n";
        error_msg += "                         "+BoundaryConditionNames::periodic()+"\n";
        error_msg += "       If your boundary condition is not one of these types, then \n";
        error_msg += "       you must specify the boundary condition type for each Variable\n";
        error_msg += "       section. Please have a look at the examples.\n";
        libmesh_error_msg(error_msg);
      }
  }

  void DefaultBCBuilder::build_axisymmetric_bcs( const GetPot& input,
                                                 MultiphysicsSystem& system,
                                                 const std::set<BoundaryID>& bc_ids,
                                                 libMesh::DofMap& dof_map,
                                                 const std::string& bc_type,
                                                 std::set<std::string>& var_sections,
                                                 std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs )
  {
    // Axisymmetric depends on the variable type.
    // So, we loop over all the Variable names and prepend the type
    // with the variable name. So, the factories for axisymmetric
    // need to be instantiated with, e.g., Velocity_axisymmetric
    for( std::set<std::string>::const_iterator vars = var_sections.begin();
         vars != var_sections.end(); ++vars )
      {
        const std::string& var_section = *vars;

        // Grab FEVariable
        const FEVariablesBase& fe_var =
          GRINS::GRINSPrivate::VariableWarehouse::get_variable(var_section);

        // We don't need to do anything for constraint variables
        if( fe_var.is_constraint_var() )
          continue;

        std::string full_bc_type = var_section+"_"+bc_type;

        // Axisymmetric is Dirichlet or Neumann, depending on the Variable
        if( this->is_dirichlet_bc_type(full_bc_type) )
          {
            this->construct_dbc_core( input,
                                      system,
                                      bc_ids,
                                      fe_var,
                                      std::string("dummy"), // No section needed for axisymmetric
                                      full_bc_type,
                                      dof_map );
          }
        else if( this->is_neumann_bc_type(full_bc_type) )
          {
            this->construct_nbc_core( input,
                                      system,
                                      bc_ids,
                                      fe_var,
                                      std::string("dummy"), // No section needed for axisymmetric
                                      full_bc_type,
                                      neumann_bcs );
          }
        else
          libmesh_error_msg("ERROR! Invalid axisymmetric type: "+full_bc_type);
      }
  }

  void DefaultBCBuilder::build_bcs_by_var_section(const GetPot& input,
                                                  MultiphysicsSystem& system,
                                                  const std::string& bc_name,
                                                  const std::set<BoundaryID>& bc_ids,
                                                  libMesh::DofMap& dof_map,
                                                  std::set<std::string>& var_sections,
                                                  const std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >& bc_id_to_subdomain_id_map,
                                                  std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs)
  {
    for( std::set<std::string>::const_iterator vars = var_sections.begin();
         vars != var_sections.end(); ++vars )
      {
        const std::string& var_section = *vars;

        // Setup the input section we're in. We use this, but
        // this is also the section where the BCFactory will try
        // to parse the values, if it needs to, so let's set up
        // once and reuse it.
        std::string input_section = std::string(BoundaryConditionNames::bc_section()+"/"+bc_name+"/"+var_section);

        // All the boundary ids have the same subdomain id (this is checked
        // earlier), so just grab the first one
        std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >::const_iterator
          subdomain_ids_iter = bc_id_to_subdomain_id_map.find(*bc_ids.begin());

        // If we're on a DistributedMesh we might need to learn about
        // subdomain ids from other processors.
        std::vector<libMesh::subdomain_id_type> subdomain_ids;
        if (subdomain_ids_iter != bc_id_to_subdomain_id_map.end())
          subdomain_ids = subdomain_ids_iter->second;
        system.comm().allgather(subdomain_ids);
        std::sort( subdomain_ids.begin(), subdomain_ids.end() );
        subdomain_ids.erase
	  (std::unique(subdomain_ids.begin(), subdomain_ids.end()),
           subdomain_ids.end() );

        // Grab the FEVariable
        const FEVariablesBase& fe_var =
          GRINSPrivate::VariableWarehouse::get_variable(var_section);

        bool var_active = this->is_var_active( fe_var, subdomain_ids );

        // If the variable is *not* active and the section is there,
        // that's an error.
        if( !var_active && input.have_section(input_section) )
          {
            std::stringstream error_msg;
            error_msg << "ERROR: Cannot specify boundary condition for variable "
                      << var_section << " on boundary " << bc_name << std::endl
                      << "since it is inactive on the subdomain associated "
                      << "with this boundary." <<  std::endl;
            libmesh_error_msg(error_msg.str());
          }

        // Make sure this section is there,
        // unless that variable is a constraint variable
        // or it's not enabled on this subdomain
        if( !input.have_section(input_section) )
          {
            if( fe_var.is_constraint_var() || !var_active )
              continue;
            else
              libmesh_error_msg("ERROR: Could not find boundary condition specification for "+input_section+"!");

          }

        // Grab the type of the boundary condition
        // There may be more than one type (e.g. pin displacement in
        // two directions and load in the third direction).
        std::string bc_type_section = input_section+"/type";
        unsigned int n_bc_types = input.vector_variable_size(bc_type_section);

        for( unsigned int bc_type_idx = 0; bc_type_idx < n_bc_types; bc_type_idx++ )
          {
            std::string bc_type = input(bc_type_section, std::string("DIE!"), bc_type_idx);

            if( this->is_dirichlet_bc_type(bc_type) )
              {
                this->construct_dbc_core( input,system, bc_ids,
                                          fe_var, input_section,
                                          bc_type, dof_map );
              }
            else if( this->is_neumann_bc_type(bc_type) )
              {
                this->construct_nbc_core( input,system, bc_ids,
                                          fe_var, input_section,
                                          bc_type, neumann_bcs );
              }
            else
              libmesh_error_msg("ERROR: Invalid bc_type "+bc_type+"!");

          } // end loop over bc_types

      } // end loop over variable sections
  }

  void DefaultBCBuilder::parse_and_build_bc_id_map( const GetPot& input,
                                                    std::map<std::string,std::set<BoundaryID> >& bc_id_map )
  {
    // First, make sure the proper input variables have been set
    if( !input.have_variable(BoundaryConditionNames::bc_ids_var()) ||
        !input.have_variable(BoundaryConditionNames::bc_id_name_map_var()) )
      {
        std::string error_msg = "ERROR: Must specify both "+BoundaryConditionNames::bc_ids_var();
        error_msg += " and "+BoundaryConditionNames::bc_id_name_map_var()+"!";
        libmesh_error_msg(error_msg);
      }

    // Make sure the vectors are the same size
    unsigned int n_ids = input.vector_variable_size(BoundaryConditionNames::bc_ids_var());
    unsigned int n_names = input.vector_variable_size(BoundaryConditionNames::bc_id_name_map_var());

    if( n_ids != n_names )
      libmesh_error_msg("ERROR: Must have matching number of boundary id sets and boundary names!");

    // We'll build this up along the way and then double check
    // that we have a one-to-one match with the boundary ids
    // that the mesh knows about
    std::set<BoundaryID> all_bc_ids;

    // Now build up the map
    for( unsigned int i = 0; i < n_names; i++ )
      {
        std::string bc_name = input(BoundaryConditionNames::bc_id_name_map_var(),std::string("DIE!"),i);
        std::string bc_ids_input = input(BoundaryConditionNames::bc_ids_var(),std::string("DIE!"),i);

        // Users can group multiple bc_ids together with the ':' delimiter
        std::vector<std::string> bc_ids_str;
        StringUtilities::split_string( bc_ids_input, ":", bc_ids_str );
        std::set<BoundaryID> bc_ids;

        for(std::vector<std::string>::const_iterator it = bc_ids_str.begin();
            it < bc_ids_str.end(); ++it )
          {
            BoundaryID id = StringUtilities::string_to_T<BoundaryID>(*it);

            // Can only specify a boundary id once
            if( (bc_ids.find(id) != bc_ids.end()) ||
                (all_bc_ids.find(id) != all_bc_ids.end())   )
              libmesh_error_msg("ERROR: Can only specify a boundary ID once!");

            bc_ids.insert(id);
            all_bc_ids.insert(id);
          }

        // Insert pair into the bc_id_map
        bc_id_map.insert( std::make_pair( bc_name, bc_ids ) );
      }
  }

  void DefaultBCBuilder::verify_bc_ids_with_mesh( const MultiphysicsSystem& system,
                                                  const std::map<std::string,std::set<BoundaryID> >& bc_id_map ) const
  {
    const libMesh::MeshBase& mesh = system.get_mesh();
    const libMesh::BoundaryInfo& boundary_info = mesh.get_boundary_info();

    std::set<BoundaryID> mesh_ids = boundary_info.get_boundary_ids();
    mesh.comm().set_union(mesh_ids);

    // Collect all the bc_ids into one set so we can just compare the sets
    std::set<BoundaryID> all_ids;

    for( std::map<std::string,std::set<BoundaryID> >::const_iterator bc_it = bc_id_map.begin();
         bc_it != bc_id_map.end(); ++bc_it )
      all_ids.insert( (bc_it->second).begin(), (bc_it->second).end() );

    if( mesh_ids != all_ids )
      {
        std::string err_msg = "ERROR: Mismatch between specified boundary ids and the boundary ids in the mesh!\n";
        err_msg += "User specified ids: ";

        for( std::set<BoundaryID>::const_iterator it = all_ids.begin();
             it != all_ids.end(); ++it )
          err_msg += StringUtilities::T_to_string<BoundaryID>(*it)+" ";

        err_msg += "\n";

        err_msg += "Mesh specified ids: ";

        for( std::set<BoundaryID>::const_iterator it = mesh_ids.begin();
             it != mesh_ids.end(); ++it )
          err_msg += StringUtilities::T_to_string<BoundaryID>(*it)+" ";

        err_msg += "\n";

        libmesh_error_msg(err_msg);
      }
  }

  void DefaultBCBuilder::build_periodic_bc( const GetPot& input,
                                            libMesh::System& system,
                                            const std::set<BoundaryID>& bc_ids,
                                            const std::string& section )
  {
    // Make sure we're not using a ParallelMesh
    // https://github.com/libMesh/libmesh/issues/977
    const libMesh::MeshBase& mesh = system.get_mesh();
    const libMesh::ParallelMesh* pmesh = dynamic_cast<const libMesh::ParallelMesh*>(&mesh);
    if(pmesh)
      {
        std::stringstream error_msg;
        error_msg << "ERROR: Cannot use ParallelMesh with periodic boundary conditions!"
                  << std::endl
                  << "       See https://github.com/libMesh/libmesh/issues/977 for discussion."
                  << std::endl;
        libmesh_error_msg(error_msg.str());
      }

    libMesh::boundary_id_type invalid_bid =
      std::numeric_limits<libMesh::boundary_id_type>::max();

    libMesh::boundary_id_type slave_id = invalid_bid;
    libMesh::boundary_id_type master_id = invalid_bid;

    this->parse_periodic_master_slave_ids(input,section,master_id,slave_id);

    if( bc_ids.find(slave_id) == bc_ids.end() ||
        bc_ids.find(master_id) == bc_ids.end() )
      libmesh_error_msg("ERROR: Mismatch between bc_ids and master/slave ids for perioid bcs!");

    libMesh::RealVectorValue offset_vector =
      this->parse_periodic_offset(input,section);

    libMesh::DofMap& dof_map = system.get_dof_map();

    this->add_periodic_bc_to_dofmap( master_id, slave_id, offset_vector, dof_map );
  }

  void DefaultBCBuilder::parse_periodic_master_slave_ids
  ( const GetPot& input, const std::string& section,
    libMesh::boundary_id_type& master_id,
    libMesh::boundary_id_type& slave_id ) const
  {
    if( !input.have_variable(section+"/master_id") )
      libmesh_error_msg("ERROR: Could not find master_id in section "+section);

    if( !input.have_variable(section+"/slave_id") )
      libmesh_error_msg("ERROR: Could not find slave_id in section "+section);

    libMesh::boundary_id_type invalid_bid =
      std::numeric_limits<libMesh::boundary_id_type>::max();

    master_id = input(section+"/master_id", invalid_bid);
    slave_id = input(section+"/slave_id", invalid_bid);
  }

  libMesh::RealVectorValue DefaultBCBuilder::parse_periodic_offset
  (const GetPot& input, const std::string& section) const
  {
    std::string input_section = section+"/boundary_offset";
    if( !input.have_variable(input_section) )
      libmesh_error_msg("ERROR: Could not find boundary_offset in section "+section);

    unsigned int n_comps = input.vector_variable_size(input_section);

    if( n_comps > 3 )
      libmesh_error_msg("ERROR: Cannot specify more than 3 components for boundary_offset!");

    libMesh::RealVectorValue offset;
    libMesh::Real invalid_real = std::numeric_limits<libMesh::Real>::max();

    for( unsigned int i = 0; i < n_comps; i++ )
      offset(i) = input(input_section, invalid_real, i );

    return offset;
  }

  void DefaultBCBuilder::build_bc_to_subdomain_map_check_with_mesh
  ( const MultiphysicsSystem& system,
    std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >& bc_id_to_subdomain_id_map ) const
  {
    std::vector<libMesh::dof_id_type> elem_ids;
    std::vector<unsigned short int> side_ids;
    std::vector<BoundaryID> bc_ids;

    const libMesh::MeshBase& mesh = system.get_mesh();

    // Extract the list of elements on the boundary
    mesh.get_boundary_info().build_active_side_list( elem_ids, side_ids, bc_ids );

    libmesh_assert_equal_to( elem_ids.size(), side_ids.size() );
    libmesh_assert_equal_to( elem_ids.size(), bc_ids.size() );

    for( unsigned int i = 0; i < elem_ids.size(); i++ )
      {
        const libMesh::Elem* elem_ptr = mesh.elem(elem_ids[i]);

        libMesh::subdomain_id_type subdomain_id = elem_ptr->subdomain_id();

        std::map<BoundaryID,std::vector<libMesh::subdomain_id_type> >::iterator
          map_it = bc_id_to_subdomain_id_map.find(bc_ids[i]);

        if( map_it == bc_id_to_subdomain_id_map.end() )
          {
            std::vector<libMesh::subdomain_id_type> sid(1,subdomain_id);
            bc_id_to_subdomain_id_map.insert(std::make_pair(bc_ids[i],sid));
          }
        else
          {
            std::vector<libMesh::subdomain_id_type>& sids = map_it->second;
            bool found_sid = false;

            for( unsigned int i = 0; i < sids.size(); i++ )
              {
                if( sids[i] == subdomain_id )
                  found_sid = true;
              }

            if( !found_sid && (sids.size() > 2) )
              libmesh_error_msg("ERROR: How do you have more than 2 subdomain ids for a boundary id?!");
            else if( !found_sid && sids.size() < 2 )
              sids.push_back(subdomain_id);
          }
      }
  }

  bool DefaultBCBuilder::is_var_active( const FEVariablesBase& var,
                                        const std::vector<libMesh::subdomain_id_type>& subdomain_ids ) const
  {
    bool var_active = false;

    // If the var subdomain_ids are empty, then this var is active
    // on the whole domain and we don't need to do anything
    if( var.subdomain_ids().empty() )
      var_active = true;

    else
      {
        // Now check if this variable is enabled on this subdomain
        const std::set<libMesh::subdomain_id_type>& var_subdomain_ids =
          var.subdomain_ids();

        for(std::vector<libMesh::subdomain_id_type>::const_iterator id = subdomain_ids.begin(); id < subdomain_ids.end(); ++id )
          {
            if( var_subdomain_ids.find(*id) != var_subdomain_ids.end() )
              var_active = true;
          }
      }

    return var_active;
  }

} // end namespace GRINS
