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
#include "grins/electrostatics.h"

// GRINS
#include "grins/electrostatics_bc_handling.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  Electrostatics::Electrostatics( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _V_var_name( input("Physics/VariableNames/ElectricPotential", GRINS::V_var_name_default ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+electrostatics+"/FE_family", "LAGRANGE") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+electrostatics+"/V_order", "SECOND") ) )
  {
    //_bc_handler = new GRINS::ElectrostaticsBCHandling( physics_name, input );

    return;
  }

  Electrostatics::~Electrostatics()
  {
    return;
  }

  void Electrostatics::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _V_var = system->add_variable( _V_var_name, _V_order, _V_FE_family);

    return;
  }

  void Electrostatics::init_context( libMesh::FEMContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.element_fe_var[_V_var]->get_JxW();
    context.element_fe_var[_V_var]->get_phi();
    context.element_fe_var[_V_var]->get_dphi();
    context.element_fe_var[_V_var]->get_xyz();
    
    context.side_fe_var[_V_var]->get_JxW();
    context.side_fe_var[_V_var]->get_phi();
    context.side_fe_var[_V_var]->get_dphi();
    context.side_fe_var[_V_var]->get_xyz();
    
    return;
  }

  void Electrostatics::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    system->time_evolving(_V_var);

    return;
  }

  void Electrostatics::element_time_derivative( bool compute_jacobian,
						libMesh::FEMContext& context )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("Electrostatics::element_time_derivative");
#endif
    
    // The number of local degrees of freedom in each variable.
    const unsigned int n_V_dofs = context.dof_indices_var[_V_var].size();
  
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_V_var]->get_JxW();
    
    // The electric potential shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi =
      context.element_fe_var[_V_var]->get_dphi();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<Number> &F_V = *context.elem_subresiduals[_V_var]; // R_{V}

    libMesh::DenseSubMatrix<Number> &K_VV = *context.elem_subjacobians[_V_var][_V_var]; // R_{V},{V}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number V;
	libMesh::Gradient grad_V;

	context.interior_value(_V_var, qp, V);
	context.interior_gradient(_V_var, qp, grad_V);
	
	// Loop over electric potential dofs
	for (unsigned int i=0; i != n_V_dofs; i++)
	  {
	    F_V(i) += grad_V*V_gradphi[i][qp]*JxW[qp];
	    
	    if( compute_jacobian )
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    K_VV(i,j) += V_gradphi[j][qp]*V_gradphi[i][qp]*JxW[qp];
		  }
	      }
	  } // electric potential dofs

      } // end quadrature loop

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("Electrostatics::element_time_derivative");
#endif

    return;
  }


  void Electrostatics::side_time_derivative( bool compute_jacobian,
					     libMesh::FEMContext& context )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("Electrostatics::side_time_derivative");
#endif

    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
         it != ids.end(); it++ )
      {
        libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);
	
	_bc_handler->apply_neumann_bcs( context, _V_var, compute_jacobian, *it );
      }

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("Electrostatics::side_time_derivative");
#endif
    
    return;
  }

} // namespace GRINS
