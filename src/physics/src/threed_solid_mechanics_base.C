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

// This class
#include "grins/threed_solid_mechanics_base.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  ThreeDSolidMechanicsBase::ThreeDSolidMechanicsBase( const PhysicsName & physics_name,
                                                      const PhysicsName & core_physics_name,
                                                      const GetPot & input )
    :  SolidMechanicsAbstract<3>(physics_name,core_physics_name,input)
  {
    if( this->_disp_vars.dim() != 3 )
      {
        std::string msg = "ERROR: "+physics_name+" only valid for three dimensions!\n";
        msg += "       Make sure you have three components in your Displacement type variable.\n";
        libmesh_error_msg(msg);
      }
  }

  void ThreeDSolidMechanicsBase::init_context( AssemblyContext & context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_dphi();
    this->get_fe(context)->get_xyz();
  }

  libMesh::Tensor ThreeDSolidMechanicsBase::form_def_gradient( const libMesh::Gradient & grad_u,
                                                               const libMesh::Gradient & grad_v,
                                                               const libMesh::Gradient & grad_w ) const
  {
    return libMesh::Tensor( 1.0+grad_u(0), grad_u(1), grad_u(2),
                            grad_v(0), 1.0+grad_v(1), grad_v(2),
                            grad_w(0), grad_w(1), 1.0+grad_w(2) );
  }

  void ThreeDSolidMechanicsBase::compute_invariants( const libMesh::Tensor & C,
                                                     libMesh::Number & I1,
                                                     libMesh::Number & I2,
                                                     libMesh::Number & I3 ) const
  {
    I1 = C.tr();

    // I2 = 0.5*( (tr(C))^2 - tr(C^2) )
    I2 = 0.5*( I1*I1 - (C*C).tr() );

    I3 = C.det();
  }
} // end namespace GRINS
