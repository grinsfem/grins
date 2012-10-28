//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id:$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "axisym_electrostatics.h"
#include "axisym_electrostatics_bc_handling.h"

namespace GRINS
{
  AxisymmetricElectrostatics::AxisymmetricElectrostatics(const std::string& physics_name, const GetPot& input )
    : Physics(physics_name)
  {
    this->read_input_options(input);
  
    _bc_handler = new GRINS::AxisymmetricElectrostaticsBCHandling( physics_name, input );
    
    return;
  }

  AxisymmetricElectrostatics::~AxisymmetricElectrostatics()
  {
    return;
  }

  void AxisymmetricElectrostatics::read_input_options( const GetPot& input )
  { 
    this->_V_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_electrostatics+"/FE_family", "LAGRANGE") );
    
    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_electrostatics+"/V_order", "SECOND") );

    this->_V_var_name = input("Physics/VariableNames/ElectricPotential", GRINS::V_var_name_default );
    
    return;
  }

  void AxisymmetricElectrostatics::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _V_var = system->add_variable( _V_var_name, _V_order, _V_FE_family);

    return;
  }

  void AxisymmetricElectrostatics::init_context( libMesh::DiffContext &context )
  {
    libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);
    
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    c.element_fe_var[_V_var]->get_JxW();
    c.element_fe_var[_V_var]->get_phi();
    c.element_fe_var[_V_var]->get_dphi();
    c.element_fe_var[_V_var]->get_xyz();
    
    c.side_fe_var[_V_var]->get_JxW();
    c.side_fe_var[_V_var]->get_phi();
    c.side_fe_var[_V_var]->get_dphi();
    c.side_fe_var[_V_var]->get_xyz();
    
    return;
  }

  void AxisymmetricElectrostatics::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    
    system->time_evolving(_V_var);

    return;
  }

  bool AxisymmetricElectrostatics::element_time_derivative( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* /*system*/ )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricElectrostatics::element_time_derivative");
#endif

    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    // The number of local degrees of freedom in each variable.
    const unsigned int n_V_dofs = c.dof_indices_var[_V_var].size();
  
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      c.element_fe_var[_V_var]->get_JxW();
    
    // The electric potential shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi =
      c.element_fe_var[_V_var]->get_dphi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint =
      c.element_fe_var[_V_var]->get_xyz();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<Number> &F_V = *c.elem_subresiduals[_V_var]; // R_{V}

    libMesh::DenseSubMatrix<Number> &K_VV = *c.elem_subjacobians[_V_var][_V_var]; // R_{V},{V}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = c.element_qrule->n_points();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Number r = qpoint[qp](0);

	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number V;
	libMesh::Gradient grad_V;

	c.interior_value<Real>(_V_var, qp, V);
	c.interior_gradient<Real>(_V_var, qp, grad_V);
	
	// Loop over electric potential dofs
	for (unsigned int i=0; i != n_V_dofs; i++)
	  {
	    F_V(i) += grad_V*V_gradphi[i][qp]*r*JxW[qp];
	    
	    if( request_jacobian )
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    K_VV(i,j) += V_gradphi[j][qp]*V_gradphi[i][qp]*r*JxW[qp];
		  }
	      }
	  } // electric potential dofs

      } // end quadrature loop

    return request_jacobian;
  }

  bool AxisymmetricElectrostatics::element_constraint( bool request_jacobian,
						       libMesh::DiffContext& /*context*/,
						       libMesh::FEMSystem* /*system*/ )
  {
    return request_jacobian;
  }

  bool AxisymmetricElectrostatics::side_time_derivative( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricElectrostatics::side_time_derivative");
#endif
    
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    const GRINS::BoundaryID boundary_id =
      system->get_mesh().boundary_info->boundary_id(c.elem, c.side);
    
    libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);
    
    _bc_handler->apply_neumann_bcs( c, _V_var, request_jacobian, boundary_id );

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("AxisymmetricElectrostatics::side_time_derivative");
#endif
    
    return request_jacobian;
  }
  
  bool AxisymmetricElectrostatics::side_constraint( bool request_jacobian,
						    libMesh::DiffContext& /*context*/,
						    libMesh::FEMSystem* /*system*/ )
  {
#ifdef USE_GRVY_TIMERS
    //this->_timer->BeginTimer("AxisymmetricElectrostatics::side_constraint");
#endif
    
    //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
#ifdef USE_GRVY_TIMERS
    //this->_timer->EndTimer("AxisymmetricElectrostatics::side_constraint");
#endif

    return request_jacobian;
  }
  

  bool AxisymmetricElectrostatics::mass_residual( bool request_jacobian,
						  libMesh::DiffContext&,
						  libMesh::FEMSystem* )
  {
    return request_jacobian;
  }

} //namespace GRINS
