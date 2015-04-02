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
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_multipointer.h"

namespace GRINS
{
  void ParameterUser::set_parameter( libMesh::Number & param_variable,
                                     const GetPot& input,
                                     const std::string & param_name,
                                     libMesh::Number param_default )
  {
    param_variable = input(param_name, param_default);

    libmesh_assert_equal_to (_my_parameters.count(param_name), 0);

    _my_parameters[param_name] = &param_variable;
  }

  void ParameterUser::register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiPointer<libMesh::Number> & param_pointer )
    const
  {
    std::map<std::string, libMesh::Number*>::const_iterator it =
      _my_parameters.find(param_name);

    if (it != _my_parameters.end())
      {
        std::cout << _my_name << " uses parameter " << param_name
                  << std::endl;
        param_pointer.push_back(it->second);
      }
  }

} // namespace GRINS
