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
#include "grins/pressure_pinning.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"

namespace GRINS
{

  PressurePinning::PressurePinning( const GetPot& input,
                                    const std::string& physics_name )
    : _pinned_elem_id(libMesh::DofObject::invalid_id)
  {
    _pin_value = input("Physics/"+physics_name+"/pin_value", 0.0 );

    unsigned int pin_loc_dim = input.vector_variable_size("Physics/"+physics_name+"/pin_location");

    // If the user is specifying a pin_location, it had better be at least 2-dimensional
    if( pin_loc_dim > 0 && pin_loc_dim < 2 )
      {
        std::cerr << "Error: pressure pin location must be at least 2 dimensional"
                  << std::endl;
        libmesh_error();
      }

    _pin_location(0) = input("Physics/"+physics_name+"/pin_location", 0.0, 0 );

    if( pin_loc_dim > 1 )
      _pin_location(1) = input("Physics/"+physics_name+"/pin_location", 0.0, 1 );

    if( pin_loc_dim == 3 )
      _pin_location(2) = input("Physics/"+physics_name+"/pin_location", 0.0, 2 );
  }

  void PressurePinning::check_pin_location( const libMesh::MeshBase& mesh )
  {
    // We need to reset to invalid_id since this may not be the first time called
    _pinned_elem_id = libMesh::DofObject::invalid_id;

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
      {
        const libMesh::Elem* elem = *el;

        if( elem->contains_point(_pin_location) )
          {
            _pinned_elem_id = elem->id();
            break;
          }
      }

    // If we found the point on one of the processors, then we need
    // to tell all the others. invalid_id is exceedingly large,
    // so if we found an element, that id should be the smallest
    mesh.comm().min( _pinned_elem_id );

    if( _pinned_elem_id == libMesh::DofObject::invalid_id )
      {
        libMesh::err << "ERROR: Could not locate point " << _pin_location
                     << " in mesh!" << std::endl;
        libmesh_error();
      }
  }

  void PressurePinning::pin_value( libMesh::DiffContext &context,
                                   const bool request_jacobian,
                                   const VariableIndex var,
                                   const double penalty )
  {
    // Make sure we've called check_pin_location() and that pin location
    // is in the mesh somewhere
    libmesh_assert_not_equal_to( _pinned_elem_id, libMesh::DofObject::invalid_id );

    /** \todo pin_location needs to be const. Currently a libMesh restriction. */
    AssemblyContext &c = libMesh::cast_ref<AssemblyContext&>(context);

    if( c.get_elem().id() == _pinned_elem_id )
      {
        // This is redundant for vast majority of cases, but we trying to
        // be prepared for cases involving, e.g. mesh motion that we're not
        // currently handling.
        if( !c.get_elem().contains_point(_pin_location) )
          {
            libmesh_error_msg("ERROR: _pin_location not in the current element!");
          }

        libMesh::DenseSubVector<libMesh::Number> &F_var = c.get_elem_residual(var); // residual
        libMesh::DenseSubMatrix<libMesh::Number> &K_var = c.get_elem_jacobian(var, var); // jacobian

        // The number of local degrees of freedom in p variable.
        const unsigned int n_var_dofs = c.get_dof_indices(var).size();

        libMesh::Number var_value = c.point_value(var, _pin_location);

        libMesh::FEType fe_type = c.get_element_fe(var)->get_fe_type();

        libMesh::Point point_loc_in_masterelem =
          libMesh::FEInterface::inverse_map(c.get_dim(), fe_type, &c.get_elem(), _pin_location);

        std::vector<libMesh::Real> phi(n_var_dofs);

        for (unsigned int i=0; i != n_var_dofs; i++)
          phi[i] = libMesh::FEInterface::shape( c.get_dim(), fe_type, &c.get_elem(), i,
                                                point_loc_in_masterelem );

        for (unsigned int i=0; i != n_var_dofs; i++)
          {
            F_var(i) += penalty*(var_value - _pin_value)*phi[i];

            /** \todo What the hell is the c.get_elem_solution_derivative() all about? */
            if (request_jacobian && c.get_elem_solution_derivative())
              {
                libmesh_assert (c.get_elem_solution_derivative() == 1.0);

                for (unsigned int j=0; j != n_var_dofs; j++)
                  K_var(i,j) += penalty*phi[i]*phi[j];

              } // End if request_jacobian
          } // End i loop
      } // End if pin_location
  }

} // namespace GRINS
