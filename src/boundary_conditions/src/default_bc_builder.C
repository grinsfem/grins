//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
    this->parse_var_sections( input, system, var_sections );

    for( std::map<std::string,std::set<BoundaryID> >::const_iterator bc_it = bc_id_map.begin();
         bc_it != bc_id_map.end(); ++bc_it )
      {
        const std::string& bc_name = bc_it->first;

        const std::set<BoundaryID>& bc_ids = bc_it->second;

        // First check for special types of boundary conditions that can be
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
                                           var_sections,neumann_bcs);
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

} // end namespace GRINS
