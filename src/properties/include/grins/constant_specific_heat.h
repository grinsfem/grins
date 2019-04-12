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


#ifndef GRINS_CONSTANT_SPECIFIC_HEAT_H
#define GRINS_CONSTANT_SPECIFIC_HEAT_H

// GRINS
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"

class GetPot;

namespace GRINS
{
  class ConstantSpecificHeat : public ParameterUser
  {
  public:

    //! Constructor with specified material
    /*! Will look in the input file for [Materials/material/SpecificHeat/value]
      for the value of specific heat. */
    ConstantSpecificHeat( const GetPot& input, const std::string& material );

    //! Deprecated constructor
    ConstantSpecificHeat( const GetPot& input );
    virtual ~ConstantSpecificHeat();

    libMesh::Real operator()() const;

    libMesh::Real operator()( const libMesh::Real T ) const;

    libMesh::Real deriv( const libMesh::Real T ) const;

  private:

    ConstantSpecificHeat();

    libMesh::Real _cp;

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real ConstantSpecificHeat::operator()() const
  {
    return _cp;
  }

  inline
  libMesh::Real ConstantSpecificHeat::operator()( const libMesh::Real /*T*/ ) const
  {
    return (*this)();
  }

  inline
  libMesh::Real ConstantSpecificHeat::deriv( const libMesh::Real /*T*/ ) const
  {
    return 0.0;
  }

}// end namespace GRINS

#endif // GRINS_CONSTANT_SPECIFIC_HEAT_H
