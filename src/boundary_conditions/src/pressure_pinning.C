//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/pressure_pinning.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{

  PressurePinning::PressurePinning( const GetPot& input,
				    const std::string& physics_name )
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
    _pin_location(1) = input("Physics/"+physics_name+"/pin_location", 0.0, 1 );

    if( pin_loc_dim == 3 ) 
      _pin_location(2) = input("Physics/"+physics_name+"/pin_location", 0.0, 2 );

    return;
  }

  PressurePinning::~PressurePinning( )
  {
    return;
  }

  void PressurePinning::pin_value( libMesh::DiffContext &context, 
				   const bool request_jacobian,
				   const VariableIndex var, 
				   const double penalty )
  {
    /** \todo pin_location needs to be const. Currently a libMesh restriction. */
    libMesh::FEMContext &c = libMesh::libmesh_cast_ref<libMesh::FEMContext&>(context);

    if (c.elem->contains_point(_pin_location))
      {
	libMesh::DenseSubVector<libMesh::Number> &F_var = *c.elem_subresiduals[var]; // residual
	libMesh::DenseSubMatrix<libMesh::Number> &K_var = *c.elem_subjacobians[var][var]; // jacobian

	// The number of local degrees of freedom in p variable.
	const unsigned int n_var_dofs = c.dof_indices_var[var].size();

	libMesh::Number var_value = c.point_value(var, _pin_location);

	libMesh::FEType fe_type = c.element_fe_var[var]->get_fe_type();
      
	libMesh::Point point_loc_in_masterelem = 
	  libMesh::FEInterface::inverse_map(c.dim, fe_type, c.elem, _pin_location);

	std::vector<libMesh::Real> phi(n_var_dofs);

	for (unsigned int i=0; i != n_var_dofs; i++)
	  phi[i] = libMesh::FEInterface::shape( c.dim, fe_type, c.elem, i, 
						point_loc_in_masterelem );
      
	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += penalty*(var_value - _pin_value)*phi[i];
	  
	    /** \todo What the hell is the c.elem_solution_derivative all about? */
	    if (request_jacobian && c.elem_solution_derivative)
	      {
		libmesh_assert (c.elem_solution_derivative == 1.0);
	      
		for (unsigned int j=0; j != n_var_dofs; j++)
		  K_var(i,j) += penalty*phi[i]*phi[j];

	      } // End if request_jacobian
	  } // End i loop
      } // End if pin_location

    return;
  }

} // namespace GRINS
