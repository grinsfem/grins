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
#include "grins/nonlinear_solver_options.h"

namespace GRINS
{
  NonlinearSolverOptions::NonlinearSolverOptions( const GetPot & input )
    : _input(input),
      _prefix("linear-nonlinear-solver")
  {}

  void NonlinearSolverOptions::numerical_jacobian_h_vars_and_vals( std::vector<std::string> & variables,
                                                                   std::vector<libMesh::Real> & values ) const
  {
    const std::string variables_option = _prefix+"/numerical_jacobian_h_variables";
    const std::string values_option = _prefix+"/numerical_jacobian_h_values";

    const unsigned int n_numerical_jacobian_h_values =
      _input.vector_variable_size(values_option);

    // Check size consistency
    if( n_numerical_jacobian_h_values != _input.vector_variable_size(variables_option) )
      {
        std::stringstream err;
        err << "Error: found " << n_numerical_jacobian_h_values
            << " numerical_jacobian_h_values" << std::endl
            << "  but "
            << _input.vector_variable_size(variables_option)
            << " numerical_jacobian_h_variables" << std::endl;

        libmesh_error_msg(err.str());
      }

    // Now populate data, if we have any
    if( n_numerical_jacobian_h_values > 0 )
      {
        variables.resize(n_numerical_jacobian_h_values);
        values.resize(n_numerical_jacobian_h_values);

        for (unsigned int i=0; i != n_numerical_jacobian_h_values; ++i)
          {
            variables[i] = _input(variables_option,"", i);
            values[i] = _input(values_option,libMesh::Real(0), i);
          }
      }
  }

} // end namespace GRINS
