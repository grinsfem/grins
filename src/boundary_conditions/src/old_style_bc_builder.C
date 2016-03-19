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
#include "grins/old_style_bc_builder.h"

// GRINS
#include "grins/common.h"
#include "grins/var_typedefs.h"
#include "grins/physics_naming.h"
#include "grins/multiphysics_sys.h"
#include "grins/dirichlet_bc_factory_abstract.h"

// libMesh
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  void OldStyleBCBuilder::build_bcs( const GetPot& input, MultiphysicsSystem& system,
                                     std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs )
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
                                          physics_name,
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
                                          physics_name,
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
                                          physics_name,
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
                                          physics_name,
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

} // end namespace GRINS
