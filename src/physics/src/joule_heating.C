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

// This class
#include "grins/joule_heating.h"

// GRINS
#include "grins_config.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  JouleHeating::JouleHeating( const std::string& physics_name,
						      const GetPot& input )
    : Physics(physics_name,input),
      _T_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_heat_transfer+"/FE_family", "LAGRANGE") ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_electrostatics+"/FE_family", "LAGRANGE") ) ),
      _T_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_heat_transfer+"/T_order", "SECOND") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_electrostatics+"/V_order", "FIRST") ) ),
      _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) ),
      _V_var_name( input("Physics/VariableNames/ElectricPotential", V_var_name_default ) ),
      _sigma( input("Physics/"+axisymmetric_magnetostatics+"/sigma", 1.0) )
  {
    return;
  }

  JouleHeating::~JouleHeating()
  {
    return;
  }

  void JouleHeating::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, _T_order, _T_FE_family);

    _V_var  = system->add_variable(_V_var_name, _V_order, _V_FE_family);

    return;
  }

  void JouleHeating::element_time_derivative( bool compute_jacobian,
							  libMesh::FEMContext& context,
							  CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("JouleHeating::element_time_derivative");
#endif
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();
    const unsigned int n_V_dofs = context.dof_indices_var[_V_var].size();

    FEGenericBase<libMesh::Real>* T_fe;
    FEGenericBase<libMesh::Real>* V_fe;

    context.get_element_fe<libMesh::Real>( _T_var, T_fe );
    context.get_element_fe<libMesh::Real>( _V_var, V_fe );

    const std::vector<std::vector<libMesh::Real> >& T_phi = T_fe->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi = V_fe->get_dphi();

    const std::vector<libMesh::Real> &JxW = T_fe->get_JxW();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint = T_fe->get_xyz();

    libMesh::DenseSubVector<Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    libMesh::DenseSubMatrix<Number> &KTV = *context.elem_subjacobians[_T_var][_V_var]; // R_{T},{V}

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Number r = qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if( this->_is_axisymmetric )
          {
            jac *= r;
          }

	libMesh::Gradient grad_V;

	context.interior_gradient(_V_var, qp, grad_V);

	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    // Q = (j \cdot j)/\sigma = \sigma \nabla V \cdot \nabla V
	    FT(i) += _sigma*(grad_V*grad_V)*T_phi[i][qp]*jac;

	    if(compute_jacobian)
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    KTV(i,j) += _sigma*2.0*(grad_V*V_gradphi[j][qp])*T_phi[i][qp]*jac;
		  }
	      }
	  }
      } // end quadrature loop

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("JouleHeating::element_time_derivative");
#endif

    return;
  }

} // end namespace GRINS
