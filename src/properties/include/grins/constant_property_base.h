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

#ifndef GRINS_CONSTANT_PROPERTY_BASE_H
#define GRINS_CONSTANT_PROPERTY_BASE_H

//GRINS
#include "grins/property_base.h"
// C++
#include <string>

namespace GRINS
{
  //! Base class for material properties that are just constants
  /*! This class contains the basic interface and functionality. Subclasses
    should only need to handle the parsing to supply the value. */
  template<typename DerivedType>
  class ConstantPropertyBase : public PropertyBase<DerivedType>,
                               public ParameterUser
  {
  public:

    //! Constructor
    /*! parsing_key is the string used by GetPot to parse the input value.
        param_name is the name registered with ParameterUser for this quantity
        as a parameter. */
    ConstantPropertyBase( const GetPot & input,
                          const std::string & parsing_key,
                          const std::string & param_name );

    ConstantPropertyBase() = delete;

    virtual ~ConstantPropertyBase() = default;

    void reset_value( libMesh::Real new_value )
    { _value = new_value; }

    // So we can make implementations private
    friend class PropertyBase<DerivedType>;

  protected:

    libMesh::Real _value;

  private:

    libMesh::Real op_context_impl(AssemblyContext & /*context*/, unsigned int /*qp*/) const
    { return _value; }

    libMesh::Real op_point_impl(const libMesh::Point & /*p*/, const libMesh::Real /*time*/)
    { return _value; }
  };

  template<typename DerivedType>
  inline
  ConstantPropertyBase<DerivedType>::ConstantPropertyBase( const GetPot & input,
                                                           const std::string & parsing_key,
                                                           const std::string & param_name )
    : PropertyBase<DerivedType>(),
    ParameterUser(param_name),
    _value(std::numeric_limits<libMesh::Real>::max())
  {
    if(!input.have_variable(parsing_key))
      libmesh_error_msg("ERROR: Could not find input for "+parsing_key+"!\n");

    this->set_parameter(_value, input, parsing_key, _value);
  }

} // end namespace GRINS

#endif // GRINS_CONSTANT_PROPERTY_BASE_H
