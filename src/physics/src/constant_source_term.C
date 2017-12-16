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
#include "grins/constant_source_term.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  ConstantSourceTerm::ConstantSourceTerm( const std::string& physics_name, const GetPot& input )
    : SourceTermBase(physics_name,input),
      _value(0.0)
  {
    if( !input.have_variable("Physics/"+physics_name+"/Function/value") )
      {
        libMesh::err << "Error: Must specify value for ConstantSourceTerm." << std::endl
                     << "       Please specify Physics/"+physics_name+"/Function/value" << std::endl;
        libmesh_error();
      }

    this->set_parameter
      (_value, input,
       "Physics/"+physics_name+"/Function/value", _value);
  }

  void ConstantSourceTerm::element_time_derivative
  ( bool /*compute_jacobian*/, AssemblyContext & context )
  {
    for( std::vector<VariableIndex>::const_iterator v_it = _vars.begin();
         v_it != _vars.end(); ++v_it )
      {
        VariableIndex var = *v_it;

        // The number of local degrees of freedom in each variable.
        const unsigned int n_dofs = context.get_dof_indices(var).size();

        // Element Jacobian * quadrature weights for interior integration.
        const std::vector<libMesh::Real> &JxW =
          context.get_element_fe(var)->get_JxW();

        // The temperature shape functions at interior quadrature points.
        const std::vector<std::vector<libMesh::Real> >& phi =
          context.get_element_fe(var)->get_phi();

        // Get residuals
        libMesh::DenseSubVector<libMesh::Number> &F_var = context.get_elem_residual(var);

        // Now we will build the element Jacobian and residual.
        // Constructing the residual requires the solution and its
        // gradient from the previous timestep.  This must be
        // calculated at each quadrature point by summing the
        // solution degree-of-freedom values by the appropriate
        // weight functions.
        unsigned int n_qpoints = context.get_element_qrule().n_points();

        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            for (unsigned int i=0; i != n_dofs; i++)
              {
                F_var(i) += (this->_value)*phi[i][qp]*JxW[qp];
              }
          }

      } // Variable loop

    return;
  }

} // end namespace GRINS
