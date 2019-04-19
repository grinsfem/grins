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
#include "grins/overlapping_fluid_solid_map.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/unsteady_solver.h"

namespace GRINS
{
  OverlappingFluidSolidMap::OverlappingFluidSolidMap( MultiphysicsSystem & system,
                                                      const libMesh::PointLocatorBase & point_locator,
                                                      const std::set<libMesh::subdomain_id_type> & solid_ids,
                                                      const std::set<libMesh::subdomain_id_type> & fluid_ids,
                                                      const DisplacementVariable & solid_disp_vars,
                                                      bool use_old_solution )
  : _use_old_solution(use_old_solution)
  {
    if( solid_ids.empty() )
      libmesh_error_msg("ERROR: Must have at least one solid subdomain id!");

    if( fluid_ids.empty() )
      libmesh_error_msg("ERROR: Must have at least one fluid subdomain id!");

    this->build_maps(system,point_locator,solid_ids,fluid_ids,solid_disp_vars);

    this->parallel_sync(system);
  }

  void OverlappingFluidSolidMap::build_maps( MultiphysicsSystem & system,
                                             const libMesh::PointLocatorBase & point_locator,
                                             const std::set<libMesh::subdomain_id_type> & solid_ids,
                                             const std::set<libMesh::subdomain_id_type> & fluid_ids,
                                             const DisplacementVariable & solid_disp_vars )
  {
    const libMesh::MeshBase & mesh = system.get_mesh();

    libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
    libMesh::UniquePtr<libMesh::FEMContext>
      fem_context( libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release()) );

    // Swap current_local_solution with old_local_nonlinear_solution
    if(_use_old_solution)
      this->swap_old_solution(system);

    if( !mesh.is_serial() )
      libmesh_error_msg("ERROR: build_maps currently only implemented for ReplicatedMesh!");

    for( const auto & solid_subdomain_id : solid_ids )
      for( const auto & solid_elem : mesh.active_local_subdomain_elements_ptr_range(solid_subdomain_id) )
        {
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
                std::map<libMesh::dof_id_type,std::vector<unsigned int> > & map =
                  _solid_to_fluid_map[solid_elem->id()];

                // The solid quadrature point that are in this overlapping fluid/solid element pair
                std::vector<unsigned int> & solid_qps = map[fluid_elem->id()];
                solid_qps.push_back(qp);
              }

              // Now add to the fluid elem->overlapping solid elems map
              {
                std::set<libMesh::dof_id_type> & solid_set = _fluid_to_solid_map[fluid_elem->id()];
                solid_set.insert(solid_elem->id());
              }
            }
        }

    // Populate _overlapping_fluid_ids
    // This is technically redundant information, but is a (premature...)
    // optimization for fetching the list of fluid element ids later.
    for( const auto & solid_it : _solid_to_fluid_map )
      {
        libMesh::dof_id_type solid_id = solid_it.first;
        const auto & fluid_map = solid_it.second;

        std::set<libMesh::dof_id_type> fluid_ids;
        for( const auto fluid_it : fluid_map )
          fluid_ids.insert(fluid_it.first);

        _overlapping_fluid_ids.insert(std::make_pair(solid_id,fluid_ids));
      }

    // Swap back current_local_solution with old_local_nonlinear_solution
    if(_use_old_solution)
      this->swap_old_solution(system);
  }


  const std::set<libMesh::dof_id_type> & OverlappingFluidSolidMap::get_overlapping_fluid_elems
  ( const libMesh::dof_id_type solid_id ) const
  {
    const auto & it = _overlapping_fluid_ids.find(solid_id);

    if( it == _overlapping_fluid_ids.end() )
      this->map_error(solid_id,"solid");

    return it->second;
  }

  const std::set<libMesh::dof_id_type> & OverlappingFluidSolidMap::get_overlapping_solid_elems
  ( const libMesh::dof_id_type fluid_id ) const
  {
    const auto & it = _fluid_to_solid_map.find(fluid_id);

    if( it == _fluid_to_solid_map.end() )
      this->map_error(fluid_id,"fluid");

    return it->second;
  }

  const std::vector<unsigned int> & OverlappingFluidSolidMap::get_solid_qps
  ( const libMesh::dof_id_type solid_id, const libMesh::dof_id_type fluid_id ) const
  {
    const auto & solid_it = _solid_to_fluid_map.find(solid_id);

    if( solid_it == _solid_to_fluid_map.end())
      this->map_error(solid_id,"solid");

    const auto & fluid_map = solid_it->second;

    const auto & fluid_it = fluid_map.find(fluid_id);

    if( fluid_it == fluid_map.end() )
      this->map_error(fluid_id,"fluid");

    return fluid_it->second;
  }

  void OverlappingFluidSolidMap::parallel_sync( MultiphysicsSystem & system )
  {
    // This should be call on every processor
    libmesh_parallel_only(system.comm());

    const libMesh::MeshBase & mesh = system.get_mesh();

    // Parallel syncing data structure
    std::map<libMesh::processor_id_type, std::vector<std::pair<libMesh::dof_id_type,libMesh::dof_id_type>>>
      ids_to_push;

    this->pack_ids_to_push( mesh, ids_to_push );

    // Extract references to data members that we'll give to the extraction functor below
    auto & overlapping_fluid_elems = this->_overlapping_fluid_ids;
    auto & fluid_to_solid_map = this->_fluid_to_solid_map;

    // Lambda functor to extract out the sent data into the local data structures.
    // We only populate _overlapping_fluid_ids and _fluid_to_solid_map since that's
    // all we'll need on the fluid processors, we won't need quadrature points.
    auto extract_ids_functor =
    [& overlapping_fluid_elems, & fluid_to_solid_map]
    (libMesh::processor_id_type, const std::vector<std::pair<libMesh::dof_id_type,libMesh::dof_id_type>> & data)
    {
      for (auto & p : data)
        {
          const libMesh::dof_id_type solid_elem_id = p.first;
          const libMesh::dof_id_type fluid_elem_id = p.second;

          auto & fluid_set = overlapping_fluid_elems[solid_elem_id];
          fluid_set.insert(fluid_elem_id);

          auto & solid_set = fluid_to_solid_map[fluid_elem_id];
          solid_set.insert(solid_elem_id);
        }
    }; // end functor

    libMesh::Parallel::push_parallel_vector_data(system.comm(), ids_to_push, extract_ids_functor);
  }

  void OverlappingFluidSolidMap::pack_ids_to_push
    ( const libMesh::MeshBase & mesh,
      std::map<libMesh::processor_id_type,
      std::vector<std::pair<libMesh::dof_id_type,libMesh::dof_id_type>>> & ids_to_push ) const
  {
    // We need to send the solid ids to all the other processors that contain the fluid elements
    // overlapped by the solid element. To do that, we need to pack the information to send.
    // The packing is mapped by processor id. So we go through each solid element, and for all
    // the fluid elements overlapped by that solid element, we figure out the processor id of that
    // fluid element and then we pack the pair of solid element id and fluid element id. We need to
    // send both so that we can keep the association on the receiving processor and not require redundant
    // searching through the PointLocator.
    //
    // FIXME: Note we are assuming that the fluid element is present on this processor so that we can grab its
    // pointer and, therefore, processor id. This will likely break on DistributedMeshes so additional
    // work will be needed to support DistributedMeshes.
    for( const auto & solid_it : _overlapping_fluid_ids )
      {
        const libMesh::dof_id_type solid_elem_id = solid_it.first;
        const auto & fluid_set = solid_it.second;

        for( const auto & fluid_elem_id : fluid_set )
          {
            // FIXME: This is where we're making an assumption about the fluid element
            //        being on this processor.
            const libMesh::Elem & fluid_elem = mesh.elem_ref(fluid_elem_id);

            const libMesh::processor_id_type fpid = fluid_elem.processor_id();

            auto & ids_pair = ids_to_push[fpid];

            ids_pair.push_back( std::make_pair(solid_elem_id,fluid_elem_id) );
          }
      }
  }

  void OverlappingFluidSolidMap::map_error(const libMesh::dof_id_type id, const std::string & type) const
  {
    std::stringstream ss;
    ss << std::string("ERROR: Could not find ")+type+std::string(" element corresponding to element id ")
       << id << std::string(" !");
    libmesh_error_msg(ss.str());
  }

  void OverlappingFluidSolidMap::swap_old_solution( MultiphysicsSystem & system )
  {
    // Extract old nonlinear solution from TimeSolver
    // Error out if this is not an UnsteadySolver
    libMesh::TimeSolver & time_solver = system.get_time_solver();

    libMesh::UnsteadySolver * unsteady_solver = dynamic_cast<libMesh::UnsteadySolver*>(&time_solver);

    if( !unsteady_solver )
      libmesh_error_msg("ERROR: Can only call swap_old_solution when using an UnsteadySolver!");

    std::swap( unsteady_solver->old_local_nonlinear_solution, system.current_local_solution );
  }

} // end namespace GRINS
