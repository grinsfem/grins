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

#ifndef GRINS_PARSED_PROPERTY_BASE_H
#define GRINS_PARSED_PROPERTY_BASE_H

//GRINS
#include "grins/assembly_context.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"

// C++
#include <string>

namespace GRINS
{
  //! Base class for material properties based on ParsedFunction
  /*! This class contains the basic interface and functionality. Subclasses
      should only need to handle the parsing of the function-string and
      create the ParsedFunction in the local _func variable. */
  class ParsedPropertyBase
  {
  public:

    ParsedPropertyBase(){};
    virtual ~ParsedPropertyBase(){};

    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

    libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time );

    virtual void init(libMesh::FEMSystem* /*system*/){};

  protected:

    //! Returns true if function string is nonzero
    /*! The ParsedFunction used is built through a string argument.
        This function checks if the string is "0". This is useful
        for cases where the function must not be zero. */
    bool check_func_nonzero( const std::string& function ) const;

    // User specified parsed function
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > _func;

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real ParsedPropertyBase::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // FIXME: We should be getting the variable index to get the qps from the context
    // not hardcode it to be 0
    const std::vector<libMesh::Point>& x = context.get_element_fe(0)->get_xyz();

    const libMesh::Point& x_qp = x[qp];

    libMesh::Number value = (*_func)(x_qp,context.time);

    return value;
  }

  inline
  libMesh::Real ParsedPropertyBase::operator()( const libMesh::Point& p, const libMesh::Real time )
  {
    return (*_func)(p,time);
  }

} // end namespace GRINS

#endif // GRINS_PARSED_PROPERTY_BASE_H
