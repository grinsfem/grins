//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "libmesh/fem_system.h"

namespace GRINS
{

  SpalartAllmarasHelper::SpalartAllmarasHelper(const GetPot& input )
    : SpalartAllmarasParameters(input),
      _flow_vars(input)
  {}

  void SpalartAllmarasHelper::init_variables( libMesh::FEMSystem* system )
  {
    this->_dim = system->get_mesh().mesh_dimension();

    this->_flow_vars.init(system);

    return;
  }

  libMesh::Real SpalartAllmarasHelper::vorticity(AssemblyContext& context, unsigned int qp) const
  {
    libMesh::Gradient grad_u, grad_v;
    grad_u = context.interior_gradient(this->_flow_vars.u_var(), qp);
    grad_v = context.interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Real vorticity_value;
    vorticity_value = fabs(grad_v(0) - grad_u(1));

    if(this->_dim == 3)
      {
        libMesh::Gradient grad_w;
        grad_w = context.interior_gradient(this->_flow_vars.w_var(), qp);

        libMesh::Real vorticity_component_0 = grad_w(1) - grad_v(2);
        libMesh::Real vorticity_component_1 = grad_u(2) - grad_v(0);

        vorticity_value += pow(pow(vorticity_component_0, 2.0) + pow(vorticity_component_1, 2.0) + pow(vorticity_value, 2.0), 0.5);
      }

    return vorticity_value;
  }

} // namespace GRINS
