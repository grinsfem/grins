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


#ifndef GRINS_OVERLAPPING_FLUID_SOLID_MAP_H
#define GRINS_OVERLAPPING_FLUID_SOLID_MAP_H

// C++
#include <map>
#include <set>
#include <vector>

// libMesh
#include "libmesh/auto_ptr.h" // UniquePtr
#include "libmesh/point_locator_base.h"
#include "libmesh/fem_context.h"

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;
  class DisplacementVariable;

  class OverlappingFluidSolidMap
  {
  public:

    OverlappingFluidSolidMap( MultiphysicsSystem & system,
                              const libMesh::PointLocatorBase & point_locator,
                              const std::set<libMesh::subdomain_id_type> & solid_ids,
                              const std::set<libMesh::subdomain_id_type> & fluid_ids,
                              const DisplacementVariable & solid_disp_vars,
                              bool use_old_solution = false );

    virtual ~OverlappingFluidSolidMap() = default;

    OverlappingFluidSolidMap() = delete;

    bool has_overlapping_fluid_elem(const libMesh::dof_id_type elem_id) const
    { return _overlapping_fluid_ids.find(elem_id) != _overlapping_fluid_ids.end(); }

    bool has_overlapping_solid_elem(const libMesh::dof_id_type elem_id) const
    { return _fluid_to_solid_map.find(elem_id) != _fluid_to_solid_map.end(); }

    //!Returns vector with the element ids of fluid elements that overlap with the given solid elem id
    const std::set<libMesh::dof_id_type> & get_overlapping_fluid_elems
    ( const libMesh::dof_id_type solid_id ) const;

    /*! Returns the quadrature point indices corresponding to those solid quadrature points
      contained in the fluid element given by fluid_id. */
    const std::vector<unsigned int> & get_solid_qps( const libMesh::dof_id_type solid_id,
                                                     const libMesh::dof_id_type fluid_id ) const;

    //!Returns vector with the element ids of solid elements that overlap with the given fluid elem id
    const std::set<libMesh::dof_id_type> & get_overlapping_solid_elems
    ( const libMesh::dof_id_type fluid_id ) const;

  private:

    void build_maps( MultiphysicsSystem & system,
                     const libMesh::PointLocatorBase & point_locator,
                     const std::set<libMesh::subdomain_id_type> & solid_ids,
                     const std::set<libMesh::subdomain_id_type> & fluid_ids,
                     const DisplacementVariable & solid_disp_vars );

    void parallel_sync( MultiphysicsSystem & system );

    void pack_ids_to_push
    ( const libMesh::MeshBase & mesh,
      std::map<libMesh::processor_id_type,
      std::vector<std::pair<libMesh::dof_id_type,libMesh::dof_id_type>>> & ids_to_push ) const;

    void map_error(const libMesh::dof_id_type id, const std::string & type) const;

    //! Vector of element ids overlapping each solid id
    std::map<libMesh::dof_id_type,std::set<libMesh::dof_id_type>> _overlapping_fluid_ids;

    std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > >
    _solid_to_fluid_map;

    std::map<libMesh::dof_id_type,std::set<libMesh::dof_id_type>> _fluid_to_solid_map;

    bool _use_old_solution;
  };

} // end namespace GRINS

#endif // GRINS_OVERLAPPING_FLUID_SOLID_MAP_H
