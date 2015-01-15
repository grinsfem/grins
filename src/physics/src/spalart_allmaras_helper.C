//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/spalart_allmaras_helper.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/inc_navier_stokes_bc_handling.h"
#include "grins/turbulence_models_macro.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"

namespace GRINS
{
  
  SpalartAllmarasHelper::SpalartAllmarasHelper(const GetPot& input )
    : _cb1(0.1355), // Define class variables
      _sigma(2./3.),
      _cb2(0.622),
      _cw1(),
      _kappa(0.41),
      _cv1(),
      _cv2(),
      _cv3(),
      _r_lin(),
      _c_w2(),
      _c_w3(),
      _flow_vars(input)
  {        
    return;
  }
  
  SpalartAllmarasHelper::~SpalartAllmarasHelper()
  {
    return;
  }
 
  void SpalartAllmarasHelper::init_variables( libMesh::FEMSystem* system )
  {    
    this->_dim = system->get_mesh().mesh_dimension();
        
    return;
  }
    
  libMesh::Real SpalartAllmarasHelper::_vorticity(AssemblyContext& context, unsigned int qp)
  {
    libMesh::Gradient grad_u, grad_v;
    grad_u = context.interior_gradient(this->_flow_vars.u_var(), qp);
    grad_v = context.interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Real _vorticity_value;
    _vorticity_value = fabs(grad_v(0) - grad_u(1));

    if(this->_dim == 3)
      {
	libMesh::Gradient grad_w;
	grad_w = context.interior_gradient(this->_flow_vars.w_var(), qp);
	
	libMesh::Real _vorticity_component_0 = grad_w(1) - grad_v(2);
	libMesh::Real _vorticity_component_1 = grad_u(2) - grad_v(0);

	_vorticity_value += pow(pow(_vorticity_component_0, 2.0) + pow(_vorticity_component_1, 2.0) + pow(_vorticity_value, 2.0), 0.5);
      }

    return _vorticity_value;
  }
  
  libMesh::Real SpalartAllmarasHelper::_source_fn(libMesh::Number nu, libMesh::Real mu, libMesh::Real wall_distance, libMesh::Real _vorticity_value)
  {
    // Step 1
    libMesh::Real _kai = nu/mu;
    
    // Step 2
    libMesh::Real _fv1 = pow(_kai, 3.0)/(pow(_kai, 3.0) + pow(this->_cv1, 3.0));

    // Step 3
    libMesh::Real _fv2 = 1 - (_kai/(1 + _kai*_fv1));

    // Step 4
    libMesh::Real _S_bar = nu/(pow(_kappa, 2.0) * pow(wall_distance, 2.0)) ;

    // Step 5, the absolute value of the vorticity
    libMesh::Real _S = _vorticity_value;

    // Step 6
    libMesh::Real _S_tilde = 0.0;
    if(_S_bar >= -this->_cv2*_S)
      {
	_S_tilde = _S + _S_bar;
      }
    else
      {
	_S_tilde = _S + (_S*(pow(this->_cv2,2.0)*_S + this->_cv3*_S_bar))/((this->_cv3 - (2*this->_cv2))*_S - _S_bar);
      }
    
    return _S_tilde;
  }
  
  libMesh::Real SpalartAllmarasHelper::_destruction_fn(libMesh::Number nu, libMesh::Real wall_distance, libMesh::Real _S_tilde)
  {
    // Step 1
    libMesh::Real _r = 0.0;
    
    if(nu/(_S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0)) < this->_r_lin)
      {
	_r = nu/(_S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0));
      }
    else
      {
	_r = this->_r_lin;
      }

    // Step 2
    libMesh::Real _g = _r + this->_c_w2*(pow(_r,6.0) - _r);

    // Step 3
    libMesh::Real _fw = _g*pow((1 + pow(this->_c_w3,6.0))/(pow(_g,6.0) + pow(this->_c_w3,6.0)), 1.0/6.0);
    
    return _fw;
  }
  

  
} // namespace GRINS

