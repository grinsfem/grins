//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "grins/linear_reaction.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/multiphysics_sys.h"
#include "grins/physics_naming.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  LinearReaction::LinearReaction( const PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _coeff(input("Physics/"+PhysicsNaming::linear_reaction()+"/coeff", libMesh::Real(1)))
  {
    // This is deleted in the base class
    this->_ic_handler = new GenericICHandler( physics_name, input );

    std::string varname_str = "Physics/"+PhysicsNaming::linear_reaction()+"/variable";
    if (input.have_variable(varname_str))
      {
        _variablename = input(varname_str, std::string());
      }
    else
      {
        std::cerr << "Error: No variable name set for "+varname_str << std::endl;
        libmesh_error();
      }
  }

  void LinearReaction::auxiliary_init( MultiphysicsSystem & system )
  {
    _u_var = system.variable_number(_variablename);
  }

  void LinearReaction::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_u_var)->get_JxW();
    context.get_element_fe(_u_var)->get_phi();

    context.get_side_fe(_u_var)->get_nothing();
  }

  void LinearReaction::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_u_var)->get_JxW();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& phi =
      context.get_element_fe(_u_var)->get_phi();

    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(_u_var, _u_var);

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_u_var);

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number u = context.interior_value(_u_var, qp);

        // First, an i-loop over the  degrees of freedom.
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            const libMesh::Number JxWxC = JxW[qp]*_coeff;
            Fu(i) += JxWxC*(phi[i][qp]*u);

            if (compute_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kuu(i,j) += JxWxC * context.get_elem_solution_derivative() *
                  (phi[i][qp]*phi[j][qp]);
          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  }
} // namespace GRINS
