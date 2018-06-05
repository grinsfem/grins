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
                              const DisplacementVariable & solid_disp_vars );

    virtual ~OverlappingFluidSolidMap(){};

    const std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > > &
    solid_map() const
    { return _solid_to_fluid_map; }

    const std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > > &
    fluid_map() const
    { return _fluid_to_solid_map; }

  private:

    OverlappingFluidSolidMap();

    void build_maps( MultiphysicsSystem & system,
                     const libMesh::PointLocatorBase & point_locator,
                     const std::set<libMesh::subdomain_id_type> & solid_ids,
                     const std::set<libMesh::subdomain_id_type> & fluid_ids,
                     const DisplacementVariable & solid_disp_vars );

    std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > >
    _solid_to_fluid_map;

    std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > >
    _fluid_to_solid_map;
  };

} // end namespace GRINS

#endif // GRINS_OVERLAPPING_FLUID_SOLID_MAP_H
