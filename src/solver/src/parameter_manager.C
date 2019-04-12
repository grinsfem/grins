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
#include "grins/parameter_manager.h"

// GRINS
#include "grins/composite_qoi.h"
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/getpot.h"
#include "libmesh/parameter_multiaccessor.h"


namespace GRINS
{
  void ParameterManager::initialize
  ( const GetPot& input,
    const std::string & parameters_varname,
    MultiphysicsSystem & system,
    CompositeQoI * qoi)
  {
    const unsigned int n_parameters =
      input.vector_variable_size(parameters_varname);
    libmesh_assert(n_parameters);

    this->parameter_name_list.resize(n_parameters);
    this->parameter_vector.clear();
    for (unsigned int i=0; i != n_parameters; ++i)
      {
        std::string param_name =
          input(parameters_varname, std::string(), i);

        this->parameter_name_list[i] = param_name;

        libMesh::ParameterMultiAccessor<libMesh::Number> *next_param =
          new libMesh::ParameterMultiAccessor<libMesh::Number>();

        // We always have Physics solving for u
        system.register_parameter(param_name, *next_param);

        // We don't always have QoIs when solving for du/dp
        if (qoi)
          qoi->register_parameter(param_name, *next_param);

        if (next_param->size() == 0)
          {
            std::cout << "No parameters named " << param_name <<
              " found in active Physics or QoIs" << std::endl;
            libmesh_error();
          }

        this->parameter_vector.push_back
          (std::unique_ptr<libMesh::ParameterAccessor<libMesh::Number> >
           (next_param));
      }
  }

} // namespace GRINS
