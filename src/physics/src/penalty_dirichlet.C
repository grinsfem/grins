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
#include "grins/penalty_dirichlet.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  using namespace libMesh;

  PenaltyDirichlet::PenaltyDirichlet( const PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _penalty(input("Physics/"+PhysicsNaming::penalty_dirichlet()+"/penalty", Real(1e10))),
      _target(input("Physics/"+PhysicsNaming::penalty_dirichlet()+"/pin_target", Real(0)))
  {
    std::string pin_variable_str = "Physics/"+PhysicsNaming::penalty_dirichlet()+"/pin_variable";
    if (input.have_variable(pin_variable_str))
      {
        _variablename_to_pin = input(pin_variable_str, std::string());
      }
    else
      {
        std::cerr << "Error: No variable name set for "+pin_variable_str << std::endl;
        libmesh_error();
      }

    this->init_bcids(input, "Physics/"+PhysicsNaming::penalty_dirichlet());
  }

  void PenaltyDirichlet::auxiliary_init( MultiphysicsSystem & system )
  {
    _variable_to_pin =
      system.variable_number(_variablename_to_pin);
  }

  void PenaltyDirichlet::init_context( AssemblyContext & context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    context.get_side_fe(_variable_to_pin)->get_JxW();
    context.get_side_fe(_variable_to_pin)->get_phi();
  }

  void PenaltyDirichlet::side_constraint
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    if (!this->is_on_active_boundary(context))
      return;

    FEBase * side_fe = nullptr;
    context.get_side_fe(_variable_to_pin, side_fe);

    const std::vector<Real> & JxW = side_fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = side_fe->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_dofs = context.n_dof_indices(_variable_to_pin);

    // The subvectors and submatrices we need to fill:
    DenseSubMatrix<Number> & K = context.get_elem_jacobian(_variable_to_pin, _variable_to_pin);
    DenseSubVector<Number> & F = context.get_elem_residual(_variable_to_pin);

    const unsigned int n_qpoints = context.get_side_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution at the old Newton iterate
        Number u = context.side_value(_variable_to_pin, qp);

        const Real JxWxP = -JxW[qp] * _penalty;

        // The residual from the boundary terms, penalize non-target // values
        for (unsigned int i=0; i != n_dofs; i++)
          F(i) += JxWxP * (u - _target) * phi[i][qp];

        if (compute_jacobian)
          for (unsigned int i=0; i != n_dofs; i++)
            for (unsigned int j=0; j != n_dofs; ++j)
              K(i,j) += JxWxP * (phi[i][qp] * phi[j][qp]);
      } // end of the quadrature point qp-loop
  }

} // namespace GRINS
