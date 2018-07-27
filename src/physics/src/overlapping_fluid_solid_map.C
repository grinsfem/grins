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
#include "grins/overlapping_fluid_solid_map.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/elem.h"

namespace GRINS
{
  OverlappingFluidSolidMap::OverlappingFluidSolidMap( MultiphysicsSystem & system,
                                                      const libMesh::PointLocatorBase & point_locator,
                                                      const std::set<libMesh::subdomain_id_type> & solid_ids,
                                                      const std::set<libMesh::subdomain_id_type> & fluid_ids,
                                                      const DisplacementVariable & solid_disp_vars )
  {
    if( solid_ids.empty() )
      libmesh_error_msg("ERROR: Must have at least one solid subdomain id!");

    if( fluid_ids.empty() )
      libmesh_error_msg("ERROR: Must have at least one fluid subdomain id!");

    this->build_maps(system,point_locator,solid_ids,fluid_ids,solid_disp_vars);
  }

  void OverlappingFluidSolidMap::build_maps( MultiphysicsSystem & system,
                                             const libMesh::PointLocatorBase & point_locator,
                                             const std::set<libMesh::subdomain_id_type> & solid_ids,
                                             const std::set<libMesh::subdomain_id_type> & fluid_ids,
                                             const DisplacementVariable & solid_disp_vars )
  {
    const libMesh::MeshBase & mesh = system.get_mesh();

    libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
    libMesh::UniquePtr<libMesh::FEMContext> fem_context( libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release()) );

    if( !mesh.is_serial() )
      libmesh_error_msg("ERROR: build_maps currently only implemented for ReplicatedMesh!");

    for( std::set<libMesh::subdomain_id_type>::const_iterator solid_id_it = solid_ids.begin();
         solid_id_it != solid_ids.end(); ++solid_id_it )
      for( libMesh::MeshBase::const_element_iterator e = mesh.active_subdomain_elements_begin(*solid_id_it);
           e != mesh.active_local_subdomain_elements_end(*solid_id_it);
           ++e )
        {
          // Convenience
          const libMesh::Elem * solid_elem = *e;

          // Setup FEMContext for computing solid displacements
          const std::vector<libMesh::Point>& qpoints =
            fem_context->get_element_fe(solid_disp_vars.u(),2)->get_xyz();

          fem_context->get_element_fe(solid_disp_vars.u(),2)->get_phi();

          fem_context->pre_fe_reinit(system,solid_elem);
          fem_context->elem_fe_reinit();

          // Find what fluid element contains each of the quadrature points and cache
          for( unsigned int qp = 0; qp < qpoints.size(); qp++ )
            {
              libMesh::Real u_disp = 0;
              libMesh::Real v_disp = 0;
              libMesh::Real w_disp = 0;

              fem_context->interior_value(solid_disp_vars.u(), qp, u_disp);
              if( solid_disp_vars.dim() >= 2 )
                fem_context->interior_value(solid_disp_vars.v(), qp, v_disp);
              if( solid_disp_vars.dim() == 3 )
                fem_context->interior_value(solid_disp_vars.w(), qp, w_disp);

              libMesh::Point U( u_disp, v_disp, w_disp );

              // We need to look for overlap with *displaced* solid point
              libMesh::Point x = qpoints[qp]+U;

              const libMesh::Elem * fluid_elem = point_locator( x, &fluid_ids );

              if( !fluid_elem )
                libmesh_error_msg("ERROR: Could not find fluid element for given displacement! Likely solid displaced off the fluid mesh!");

              // Now add to the solid elem->overlapping fluid elems map, but only if
              // the solid elem belongs to this processor
              {
                std::map<libMesh::dof_id_type,std::vector<unsigned int> > & fluid_elem_map =
                  _solid_to_fluid_map[solid_elem->id()];

                // The solid quadrature point that are in this overlapping fluid/solid element pair
                std::vector<unsigned int>& solid_qps = fluid_elem_map[fluid_elem->id()];
                solid_qps.push_back(qp);
              }

              // Now add to the fluid elem->overlapping solid elems map, but only if
              // the fluid elem belongs to this processor
              {
                std::map<libMesh::dof_id_type,std::vector<unsigned int> > & solid_elem_map =
                  _fluid_to_solid_map[fluid_elem->id()];

                // The solid quadrature point that are in this overlapping fluid/solid element pair
                std::vector<unsigned int>& solid_qps = solid_elem_map[solid_elem->id()];
                solid_qps.push_back(qp);
              }
            }
        }
  }

} // end namespace GRINS
