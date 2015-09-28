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

#ifndef GRINS_PARAMETER_USER_H
#define GRINS_PARAMETER_USER_H

// C++
#include <string>
#include <set>

//GRINS
#include "grins_config.h"

//libMesh
#include "libmesh/libmesh.h"

//C++
#include <map>
#include <string>

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  template <typename Scalar>
  class ParameterMultiAccessor;
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
    //file and prepare it for registration with this call.
    virtual void set_parameter
      ( libMesh::Number & param_variable,
        const GetPot & input,
        const std::string & param_name,
        libMesh::Number param_default );

    // FIXME: add set_parameter for vectors

    //! Each subclass will register its copy of an independent
    //  variable when the library makes this call.
    //  If the subclass has more than one copy to register, or if the
    //  subclass needs to register a parameter which was not
    //  previously assigned with set_parameter, then this method will
    //  need to be overridden.
    virtual void register_parameter
      ( const std::string & param_name,
        libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const;

  private:
    std::map<std::string, libMesh::Number*> _my_parameters;

    // This could be more efficient as a reference now, but we'd
    // probably inadvertently break it later.
    std::string _my_name;

  }; // End Physics class declarations

  /* ------------------------- Inline Functions -------------------------*/

} // End namespace GRINS

#endif //GRINS_PARAMETER_USER_H
