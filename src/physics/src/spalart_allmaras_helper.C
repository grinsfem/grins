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
#include "grins/spalart_allmaras_helper.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/turbulence_models_macro.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/physics_naming.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  SpalartAllmarasHelper::SpalartAllmarasHelper(const GetPot& input )
    : _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,PhysicsNaming::spalart_allmaras(),VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,PhysicsNaming::spalart_allmaras(),VariablesParsing::PHYSICS)))
  {}

  libMesh::Real SpalartAllmarasHelper::vorticity(AssemblyContext& context, unsigned int qp) const
  {
    libMesh::Gradient grad_u, grad_v;
    grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
    grad_v = context.interior_gradient(this->_flow_vars.v(), qp);

    libMesh::Real vorticity_value;
    vorticity_value = fabs(grad_v(0) - grad_u(1));

    if(context.get_system().get_mesh().mesh_dimension() == 3)
      {
        libMesh::Gradient grad_w;
        grad_w = context.interior_gradient(this->_flow_vars.w(), qp);

        libMesh::Real vorticity_component_0 = grad_w(1) - grad_v(2);
        libMesh::Real vorticity_component_1 = grad_u(2) - grad_v(0);

        libMesh::Real term = vorticity_component_0*vorticity_component_0
          + vorticity_component_1*vorticity_component_1
          + vorticity_value*vorticity_value;

        vorticity_value += std::sqrt(term);
      }

    return vorticity_value;
  }

} // namespace GRINS
