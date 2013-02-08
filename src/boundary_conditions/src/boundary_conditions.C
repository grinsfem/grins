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
#include "grins/boundary_conditions.h"

// GRINS
#include "grins/cached_values.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{

  BoundaryConditions::BoundaryConditions( )
  {
    return;
  }

  BoundaryConditions::~BoundaryConditions( )
  {
    return;
  }

  void BoundaryConditions::apply_neumann( libMesh::FEMContext& context,
					  const VariableIndex var,
					  const libMesh::Real sign,
					  const libMesh::Point& value ) const
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    const std::vector<libMesh::Point> &normals = context.side_fe_var[var]->get_normals();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual

    unsigned int n_qpoints = context.side_qrule->n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*JxW_side[qp]*value*normals[qp]*var_phi_side[i][qp];
	  }
      }

    return;
  }

  void GRINS::BoundaryConditions::apply_neumann_normal( libMesh::FEMContext& context,
							const GRINS::VariableIndex var,
							const libMesh::Real sign,
							libMesh::Real value ) const
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual

    unsigned int n_qpoints = context.side_qrule->n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*value*var_phi_side[i][qp]*JxW_side[qp];
	  }
      }

    return;
  }


  void BoundaryConditions::apply_neumann_axisymmetric( libMesh::FEMContext& context,
						       const VariableIndex var,
						       const libMesh::Real sign,
						       const libMesh::Point& value ) const
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& var_qpoint =
      context.side_fe_var[var]->get_xyz();

    const std::vector<libMesh::Point> &normals = context.side_fe_var[var]->get_normals();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual

    unsigned int n_qpoints = context.side_qrule->n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Number r = var_qpoint[qp](0);

	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*r*JxW_side[qp]*value*normals[qp]*var_phi_side[i][qp];
	  }
      }

    return;
  }

  void BoundaryConditions::apply_neumann( libMesh::FEMContext& context,
					  const CachedValues& cache,
					  const bool request_jacobian,
					  const VariableIndex var,
					  const libMesh::Real sign,
					  const std::tr1::shared_ptr<NeumannFuncObj> neumann_func ) const
  {
    // The number of local degrees of freedom
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();
  
    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    const std::vector<libMesh::Point> &normals = context.side_fe_var[var]->get_normals();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual
    libMesh::DenseSubMatrix<libMesh::Number> &K_var = *context.elem_subjacobians[var][var]; // jacobian

    unsigned int n_qpoints = context.side_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Point bc_value = neumann_func->value( context, cache, qp );
	const libMesh::Point jac_value = neumann_func->derivative( context, cache, qp );

	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*JxW_side[qp]*bc_value*normals[qp]*var_phi_side[i][qp];

	    if (request_jacobian)
	      {
		for (unsigned int j=0; j != n_var_dofs; j++)
		  {
		    K_var(i,j) += sign*JxW_side[qp]*jac_value*normals[qp]*
		      var_phi_side[i][qp]*var_phi_side[j][qp];
		  }
	      }
	  }
      } // End quadrature loop

    // Now must take care of the case that the boundary condition depends on variables
    // other than var.
    std::vector<VariableIndex> other_jac_vars = neumann_func->get_other_jac_vars();

    if( (other_jac_vars.size() > 0) && request_jacobian )
      {
	for( std::vector<VariableIndex>::const_iterator var2 = other_jac_vars.begin();
	     var2 != other_jac_vars.end();
	     var2++ )
	  {
	    const unsigned int n_var2_dofs = context.dof_indices_var[*var2].size();
	    const std::vector<std::vector<libMesh::Real> >& var2_phi_side =
	      context.side_fe_var[*var2]->get_phi();

	    for (unsigned int qp=0; qp != n_qpoints; qp++)
	      {
		const libMesh::Point jac_value = neumann_func->derivative( context, cache, qp, *var2 );

		for (unsigned int i=0; i != n_var_dofs; i++)
		  {
		    for (unsigned int j=0; j != n_var2_dofs; j++)
		      {
			K_var(i,j) += sign*JxW_side[qp]*jac_value*normals[qp]*
			  var_phi_side[i][qp]*var2_phi_side[j][qp];
		      }
		  }
	      }
	  } // End loop over auxillary Jacobian variables
      }
    return;
  }

  void GRINS::BoundaryConditions::apply_neumann_normal( libMesh::FEMContext& context,
							const CachedValues& cache,
							const bool request_jacobian,
							const GRINS::VariableIndex var,
							const libMesh::Real sign,
							const std::tr1::shared_ptr<GRINS::NeumannFuncObj> neumann_func ) const
  {
    // The number of local degrees of freedom
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();
  
    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual
    libMesh::DenseSubMatrix<libMesh::Number> &K_var = *context.elem_subjacobians[var][var]; // jacobian

    unsigned int n_qpoints = context.side_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Real bc_value = neumann_func->normal_value( context, cache, qp );
	const libMesh::Real jac_value = neumann_func->normal_derivative( context, cache, qp );

	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*bc_value*var_phi_side[i][qp]*JxW_side[qp];

	    if (request_jacobian)
	      {
		for (unsigned int j=0; j != n_var_dofs; j++)
		  {
		    K_var(i,j) += sign*jac_value*
		      var_phi_side[i][qp]*var_phi_side[j][qp]*JxW_side[qp];
		  }
	      }
	  }
      } // End quadrature loop

    // Now must take care of the case that the boundary condition depends on variables
    // other than var.
    std::vector<GRINS::VariableIndex> other_jac_vars = neumann_func->get_other_jac_vars();

    if( (other_jac_vars.size() > 0) && request_jacobian )
      {
	for( std::vector<GRINS::VariableIndex>::const_iterator var2 = other_jac_vars.begin();
	     var2 != other_jac_vars.end();
	     var2++ )
	  {
	    const unsigned int n_var2_dofs = context.dof_indices_var[*var2].size();
	    const std::vector<std::vector<libMesh::Real> >& var2_phi_side =
	      context.side_fe_var[*var2]->get_phi();

	    for (unsigned int qp=0; qp != n_qpoints; qp++)
	      {
		const libMesh::Real jac_value = neumann_func->normal_derivative( context, cache, qp, *var2 );

		for (unsigned int i=0; i != n_var_dofs; i++)
		  {
		    for (unsigned int j=0; j != n_var2_dofs; j++)
		      {
			K_var(i,j) += sign*jac_value*
			  var_phi_side[i][qp]*var2_phi_side[j][qp]*JxW_side[qp];
		      }
		  }
	      }
	  } // End loop over auxillary Jacobian variables
      }
    return;
  }



  void GRINS::BoundaryConditions::apply_neumann_axisymmetric( libMesh::FEMContext& context,
							      const CachedValues& cache,
							      const bool request_jacobian,
							      const GRINS::VariableIndex var,
							      const libMesh::Real sign,
							      std::tr1::shared_ptr<GRINS::NeumannFuncObj> neumann_func ) const
  {
    // The number of local degrees of freedom
    const unsigned int n_var_dofs = context.dof_indices_var[var].size();
  
    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = context.side_fe_var[var]->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<libMesh::Real> >& var_phi_side =
      context.side_fe_var[var]->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& var_qpoint =
      context.side_fe_var[var]->get_xyz();

    const std::vector<libMesh::Point> &normals = context.side_fe_var[var]->get_normals();

    libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual
    libMesh::DenseSubMatrix<libMesh::Number> &K_var = *context.elem_subjacobians[var][var]; // jacobian

    unsigned int n_qpoints = context.side_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Point bc_value = neumann_func->value( context, cache, qp );
	const libMesh::Point jac_value = neumann_func->derivative( context, cache, qp );

	const libMesh::Number r = var_qpoint[qp](0);

	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += sign*r*JxW_side[qp]*bc_value*normals[qp]*var_phi_side[i][qp];

	    if (request_jacobian)
	      {
		for (unsigned int j=0; j != n_var_dofs; j++)
		  {
		    K_var(i,j) += sign*r*JxW_side[qp]*jac_value*normals[qp]*
		      var_phi_side[i][qp]*var_phi_side[j][qp];
		  }
	      }
	  }
      } // End quadrature loop

    // Now must take care of the case that the boundary condition depends on variables
    // other than var.
    std::vector<VariableIndex> other_jac_vars = neumann_func->get_other_jac_vars();

    if( (other_jac_vars.size() > 0) && request_jacobian )
      {
	for( std::vector<VariableIndex>::const_iterator var2 = other_jac_vars.begin();
	     var2 != other_jac_vars.end();
	     var2++ )
	  {
	    const unsigned int n_var2_dofs = context.dof_indices_var[*var2].size();
	    const std::vector<std::vector<libMesh::Real> >& var2_phi_side =
	      context.side_fe_var[*var2]->get_phi();

	    for (unsigned int qp=0; qp != n_qpoints; qp++)
	      {
		const libMesh::Number r = var_qpoint[qp](0);

		const libMesh::Point jac_value = neumann_func->derivative( context, cache, qp, *var2 );

		for (unsigned int i=0; i != n_var_dofs; i++)
		  {
		    for (unsigned int j=0; j != n_var2_dofs; j++)
		      {
			K_var(i,j) += sign*r*JxW_side[qp]*jac_value*normals[qp]*
			  var_phi_side[i][qp]*var2_phi_side[j][qp];
		      }
		  }
	      }
	  } // End loop over auxillary Jacobian variables
      }
    return;
  }

  void BoundaryConditions::pin_value( libMesh::FEMContext& context,
				      const CachedValues& /*cache*/,
				      const bool request_jacobian,
				      const VariableIndex var, 
				      const double pin_value,
				      const libMesh::Point& pin_location, 
				      const double penalty )
  {
    if (context.elem->contains_point(pin_location))
      {
	libMesh::DenseSubVector<libMesh::Number> &F_var = *context.elem_subresiduals[var]; // residual
	libMesh::DenseSubMatrix<libMesh::Number> &K_var = *context.elem_subjacobians[var][var]; // jacobian

	// The number of local degrees of freedom in p variable.
	const unsigned int n_var_dofs = context.dof_indices_var[var].size();

	libMesh::Number var_value = context.point_value(var, pin_location);

	libMesh::FEType fe_type = context.element_fe_var[var]->get_fe_type();
      
	libMesh::Point point_loc_in_masterelem = 
	  libMesh::FEInterface::inverse_map(context.dim, fe_type, context.elem, pin_location);

	std::vector<libMesh::Real> phi(n_var_dofs);

	for (unsigned int i=0; i != n_var_dofs; i++)
	  phi[i] = libMesh::FEInterface::shape( context.dim, fe_type, context.elem, i, 
						point_loc_in_masterelem );
      
	for (unsigned int i=0; i != n_var_dofs; i++)
	  {
	    F_var(i) += penalty*(var_value - pin_value)*phi[i];
	  
	    /** \todo What the hell is the context.elem_solution_derivative all about? */
	    if (request_jacobian && context.elem_solution_derivative)
	      {
		libmesh_assert (context.elem_solution_derivative == 1.0);
	      
		for (unsigned int j=0; j != n_var_dofs; j++)
		  K_var(i,j) += penalty*phi[i]*phi[j];

	      } // End if request_jacobian
	  } // End i loop
      } // End if pin_location

    return;
  }

} // namespace GRINS
