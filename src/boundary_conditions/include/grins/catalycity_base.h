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

#ifndef GRINS_CATALYCITY_BASE_H
#define GRINS_CATALYCITY_BASE_H

// C++
#include <vector>

// GRINS
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"

namespace GRINS
{
  class CatalycityBase : public ParameterUser
  {
  public:

    using ParameterUser::ParameterUser;

    virtual ~CatalycityBase() = default;

    virtual libMesh::Real operator()( const libMesh::Real T ) const =0;

    virtual libMesh::Real dT( const libMesh::Real T ) const =0;

    virtual void get_params( std::vector<libMesh::Real> & params ) =0;

    virtual void set_params( const std::vector<libMesh::Real>& params ) =0;

    //! Creates a new copy of the current class.
    /*! A raw pointer is returned and it is assumed the user will take ownership
      and worry about memory management. */
    virtual CatalycityBase* clone() const = 0;

    //! Sets parameters for use in sensitivity analysis
    virtual void set_parameters(const GetPot & input, const std::string & param_base) =0;

  };

} // end namespace GRINS

#endif // GRINS_CATALYCITY_BASE_H
