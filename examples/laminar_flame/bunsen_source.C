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
#include "bunsen_source.h"

// GRINS
#include "grins_config.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace Bunsen
{
  BunsenSource::BunsenSource( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : GRINS::Physics(physics_name,input),
      _value( input( "Physics/BunsenSource/value", 0.0 ) ),
      _r_max( input( "Physics/BunsenSource/r_max", 0.003 ) ),
      _z_min( input( "Physics/BunsenSource/z_min", 0.005 ) ),
      _z_max( input( "Physics/BunsenSource/z_max", 0.02 ) )
  {
    this->read_input_options(input);

    return;
  }

  BunsenSource::~BunsenSource()
  {
    return;
  }

  void BunsenSource::read_input_options( const GetPot& input )
  {
    this->_T_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+GRINS::reacting_low_mach_navier_stokes+"/T_FE_family", "LAGRANGE") ); 

    this->_T_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+GRINS::reacting_low_mach_navier_stokes+"/T_order", "SECOND") );

    this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

    return;
  }

  void BunsenSource::init_variables( libMesh::FEMSystem* system )
  {
    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);
    return;
  }

  void BunsenSource::element_time_derivative( bool /*compute_jacobian*/,
					      libMesh::FEMContext& context,
					      GRINS::CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("BunsenSource::element_time_derivative");
#endif
  
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_T_var]->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.element_fe_var[_T_var]->get_phi();

    // Locations of quadrature points
    const std::vector<libMesh::Point>& x_qp = context.element_fe_var[_T_var]->get_xyz();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();
    
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Real r = x_qp[qp](0);
	const libMesh::Real z = x_qp[qp](1);
	
	if( r <= _r_max &&
	    z >= _z_min &&
	    z <= _z_max )
	  {
	    for (unsigned int i=0; i != n_T_dofs; i++)
	      {
		FT(i) += _value*T_phi[i][qp]*r*JxW[qp];
	      }
	  }
      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("BunsenSource::element_time_derivative");
#endif

    return;
  }

} // namespace Bunsen
