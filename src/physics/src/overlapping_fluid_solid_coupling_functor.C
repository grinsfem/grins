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
#include "grins/overlapping_fluid_solid_coupling_functor.h"


// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

namespace GRINS
{
  void OverlappingFluidSolidCouplingFunctor::operator()
    ( const libMesh::MeshBase::const_element_iterator & range_begin,
      const libMesh::MeshBase::const_element_iterator & range_end,
      libMesh::processor_id_type p,
      map_type & coupled_elements )
  {
    for( const auto & elem : libMesh::as_range(range_begin,range_end) )
      {
        libMesh::dof_id_type elem_id = elem->id();

        bool is_solid_elem_with_overlapping_fluid_elem =
          _overlapping_map.has_overlapping_fluid_elem(elem_id);

        bool is_fluid_elem_with_overlapping_solid_elem =
          _overlapping_map.has_overlapping_solid_elem(elem_id);


        // If this element is a solid element, then we need to populate
        // the coupled_elements with all the fluid elements associated with
        // this solid element. We use the same coupling matrix for all of them.
        if( is_solid_elem_with_overlapping_fluid_elem )
          {
            const std::set<libMesh::dof_id_type> & fluid_set =
              _overlapping_map.get_overlapping_fluid_elems(elem_id);

            for( const auto & fluid_id : fluid_set )
              {
                const libMesh::Elem * fluid_elem = _mesh.elem_ptr(fluid_id);

                if(!fluid_elem)
                  libmesh_error_msg("ERROR: fluid_elem is NULL!");

                if( fluid_elem->processor_id() != p )
                  coupled_elements.insert( std::make_pair(fluid_elem,&_coupling_matrix) );
              }
          }

        // If this element is a fluid element with an overlapping solid element,
        // then we need to populate the coupled_elements with all the solid elements
        // associated with this fluid element. While we don't need the algebraic
        // coupling this will generate, it is needed to get the correct sparsity
        // pattern.
        if( is_fluid_elem_with_overlapping_solid_elem )
          {
            const std::set<libMesh::dof_id_type> & solid_set =
              _overlapping_map.get_overlapping_solid_elems(elem_id);

            for( const auto & solid_id : solid_set )
              {
                const libMesh::Elem * solid_elem = _mesh.elem_ptr(solid_id);

                if(!solid_elem)
                  libmesh_error_msg("ERROR: fluid_elem is NULL!");

                if( solid_elem->processor_id() != p )
                  coupled_elements.insert( std::make_pair(solid_elem,&_coupling_matrix) );
              }
          }

      } // end element loop
  }

} // end namespace GRINS
