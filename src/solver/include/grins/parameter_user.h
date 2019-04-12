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

#ifndef GRINS_PARAMETER_USER_H
#define GRINS_PARAMETER_USER_H

// C++
#include <string>
#include <set>

//GRINS
#include "grins_config.h"

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h" // forward declare Gradient instead?

//C++
#include <map>
#include <string>

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  template <typename Scalar>
  class ParameterMultiAccessor;

  template <typename Output, typename OutputGradient>
  class ParsedFunction;

  template <typename Output>
  class ParsedFEMFunction;
}

//! GRINS namespace
namespace GRINS
{
  //! ParameterUser base class. Utility methods for subclasses.
  /*!
    This base class defines utility methods with which a subclass can
    simultaneously set a parameter from an input file and register
    that parameter (by input file variable name) for later parameter
    sensitivity studies.
  */

  class ParameterUser
  {

  public:

    ParameterUser(const std::string & user_name) : _my_name(user_name) {}
    virtual ~ParameterUser() {}

    //! Each subclass can simultaneously read a parameter value from
    // file and prepare it for registration with this call.
    virtual void set_parameter
    ( libMesh::Number & param_variable,
      const GetPot & input,
      const std::string & param_name,
      libMesh::Number param_default );

    //! Each subclass can simultaneously read a parsed function from
    // file and prepare its inline variables for registration with this call.
    //
    // Pass a default value of "DIE!" to assert that the input file
    // contains a value for this parameter.
    virtual void set_parameter
    ( libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & func,
      const GetPot & input,
      const std::string & func_param_name,
      const std::string & param_default);

    //! Each subclass can simultaneously read a parsed function from
    // file and prepare its inline variables for registration with this call.
    //
    // Pass a default value of "DIE!" to assert that the input file
    // contains a value for this parameter.
    virtual void set_parameter
    ( libMesh::ParsedFEMFunction<libMesh::Number> & func,
      const GetPot & input,
      const std::string & func_param_name,
      const std::string & param_default);

    //! When cloning an object, we need to update parameter pointers
    //to point to the clone
    virtual void move_parameter
    (const libMesh::Number & old_parameter,
     libMesh::Number & new_parameter);

    //! When cloning an object, we need to update parameter pointers
    //to point to the clone
    virtual void move_parameter
    (const libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & old_func,
     libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient> & new_func);

    //! When cloning an object, we need to update parameter pointers
    //to point to the clone
    virtual void move_parameter
    (const libMesh::ParsedFEMFunction<libMesh::Number> & old_func,
     libMesh::ParsedFEMFunction<libMesh::Number> & new_func);

    //! A parseable function string with LIBMESH_DIM components, all 0
    static std::string zero_vector_function;

    // FIXME: add set_parameter overloads for vectors

    //! Each subclass will register its copy of an independent
    //  variable when the library makes this call.
    //  If the subclass needs to register a parameter which was not
    //  previously assigned with set_parameter, then this method will
    //  need to be overridden.
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  private:
    std::map<std::string, libMesh::Number*> _my_parameters;

    std::map<std::string,
             libMesh::ParsedFunction<libMesh::Number,libMesh::Gradient>*>
    _my_parsed_functions;

    std::map<std::string,
             libMesh::ParsedFEMFunction<libMesh::Number>*>
    _my_parsed_fem_functions;

    // This could be more efficient as a reference now, but we'd
    // probably inadvertently break it later.
    std::string _my_name;

  }; // End Physics class declarations

  /* ------------------------- Inline Functions -------------------------*/

} // End namespace GRINS

#endif //GRINS_PARAMETER_USER_H
