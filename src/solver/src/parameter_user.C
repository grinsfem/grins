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
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_multiaccessor.h"
#include "libmesh/parameter_pointer.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_fem_function_parameter.h"
#include "libmesh/parsed_function.h"
#include "libmesh/parsed_function_parameter.h"


namespace GRINS
{
#if LIBMESH_DIM == 3
  std::string ParameterUser::zero_vector_function = std::string("{0}{0}{0}");
#elif LIBMESH_DIM == 2
  std::string ParameterUser::zero_vector_function = std::string("{0}{0}");
#else
  std::string ParameterUser::zero_vector_function = std::string("{0}");
#endif


  void ParameterUser::set_parameter( libMesh::Number & param_variable,
                                     const GetPot& input,
                                     const std::string & param_name,
                                     libMesh::Number param_default )
  {
    param_variable = input(param_name, param_default);

    libmesh_assert_msg(!_my_parameters.count(param_name),
                       "ERROR: " << _my_name << " double-registered parameter " <<
                       param_name);

    _my_parameters[param_name] = &param_variable;
  }

  void ParameterUser::set_parameter
  ( libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & func,
    const GetPot & input,
    const std::string & func_param_name,
    const std::string & param_default)
  {
    if((param_default == "DIE!") &&
       (!input.have_variable(func_param_name)))
      {
        libMesh::err << "Error: Must specify parsed function for " <<
          _my_name << std::endl
                     << "       Please specify " << func_param_name << std::endl;
        libmesh_error();
      }

    func.reparse(input(func_param_name, param_default));

    libmesh_assert_msg(!_my_parsed_functions.count(func_param_name),
                       "ERROR: " << _my_name << " double-registered parameter " <<
                       func_param_name);

    _my_parsed_functions[func_param_name] = &func;
  }

  void ParameterUser::set_parameter
  ( libMesh::ParsedFEMFunction<libMesh::Number> & func,
    const GetPot & input,
    const std::string & func_param_name,
    const std::string & param_default)
  {
    if((param_default == "DIE!") &&
       (!input.have_variable(func_param_name)))
      {
        libMesh::err << "Error: Must specify parsed (fem) function for " <<
          _my_name << std::endl
                     << "       Please specify " << func_param_name << std::endl;
        libmesh_error();
      }

    func.reparse(input(func_param_name, param_default));

    libmesh_assert_msg(!_my_parsed_fem_functions.count(func_param_name),
                       "ERROR: " << _my_name << " double-registered parameter " <<
                       func_param_name);

    _my_parsed_fem_functions[func_param_name] = &func;
  }

  void ParameterUser::move_parameter
  (const libMesh::Number & old_parameter,
   libMesh::Number & new_parameter)
  {
    std::map<std::string, libMesh::Number*>::iterator it =
      _my_parameters.begin();
    const std::map<std::string, libMesh::Number*>::iterator end =
      _my_parameters.end();
    for (; it != end; ++it)
      if (it->second == &old_parameter)
        {
          it->second = &new_parameter;
          break;
        }
  }

  void ParameterUser::move_parameter
  (const libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & old_func,
   libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & new_func)
  {
    std::map
      <std::string,
       libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> *
       >::iterator it = _my_parsed_functions.begin();
    const std::map
      <std::string,
       libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> *
       >::iterator end = _my_parsed_functions.end();
    for (; it != end; ++it)
      if (it->second == &old_func)
        {
          it->second = &new_func;
          break;
        }
  }

  void ParameterUser::move_parameter
  (const libMesh::ParsedFEMFunction<libMesh::Number> & old_func,
   libMesh::ParsedFEMFunction<libMesh::Number> & new_func)
  {
    std::map
      <std::string,
       libMesh::ParsedFEMFunction<libMesh::Number> *
       >::iterator it = _my_parsed_fem_functions.begin();
    const std::map
      <std::string,
       libMesh::ParsedFEMFunction<libMesh::Number> *
       >::iterator end = _my_parsed_fem_functions.end();
    for (; it != end; ++it)
      if (it->second == &old_func)
        {
          it->second = &new_func;
          break;
        }
  }


  void ParameterUser::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    std::map<std::string, libMesh::Number*>::const_iterator it =
      _my_parameters.find(param_name);

    // Make sure we don't find duplicate parameters - a Number
    // parameter shouldn't have the same name as a ParsedFunction or
    // ParsedFEMFunction parameter.
#ifndef NDEBUG
    bool found_parameter = false;
#endif

    // First search for simple Number parameters
    if (it != _my_parameters.end())
      {
        std::cout << _my_name << " has Number parameter " << param_name
                  << std::endl;
        param_pointer.push_back
          (libMesh::ParameterPointer<libMesh::Number>(it->second));

#ifndef NDEBUG
        found_parameter = true;
#else
        return;
#endif
      }

    // Next search for inline variable parameters in parsed functions
    std::size_t last_slash_i = param_name.rfind('/');

    if (last_slash_i != std::string::npos)
      {
        std::string search_name = param_name.substr(0, last_slash_i);

        std::string var_name = param_name.substr(last_slash_i+1);

        std::map
          <std::string, libMesh::ParsedFunction
           <libMesh::Number,libMesh::Gradient>*>::const_iterator
          pf_it = _my_parsed_functions.find(search_name);

        if (pf_it != _my_parsed_functions.end())
          {
            std::cout << _my_name << " has ParsedFunction for " <<
              search_name << " / " << var_name << std::endl;
            param_pointer.push_back
              (libMesh::ParsedFunctionParameter<libMesh::Number>
               (*pf_it->second, var_name));

#ifndef NDEBUG
            libmesh_assert(!found_parameter);
            found_parameter = true;
#else
            return;
#endif
          }

        std::map
          <std::string, libMesh::ParsedFEMFunction
           <libMesh::Number>*>::const_iterator
          pff_it = _my_parsed_fem_functions.find(search_name);

        if (pff_it != _my_parsed_fem_functions.end())
          {
            std::cout << _my_name << " has ParsedFEMFunction for " <<
              search_name << " / " << var_name << std::endl;
            param_pointer.push_back
              (libMesh::ParsedFEMFunctionParameter<libMesh::Number>
               (*pff_it->second, var_name));

#ifndef NDEBUG
            libmesh_assert(!found_parameter);
#endif
          }
      }
  }

} // namespace GRINS
