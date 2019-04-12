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


#ifndef GRINS_CONSTANT_VISCOSITY_H
#define GRINS_CONSTANT_VISCOSITY_H

//GRINS
#include "grins/viscosity_base.h"
#include "grins/assembly_context.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"

#include "libmesh/fem_system.h"

class GetPot;

namespace GRINS
{
  class ConstantViscosity : public ParameterUser,
                            public ViscosityBase
  {
  public:

    //! Constructor with specified material
    /*! Will look in the input file for [Materials/material/Viscosity/value]
      for the value of viscosity. */
    ConstantViscosity( const GetPot& input, const std::string& material );

    //! Deprecated constructor
    ConstantViscosity( const GetPot& input );
    virtual ~ConstantViscosity();

    libMesh::Real operator()() const;

    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

    libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time );

    libMesh::Real operator()( const libMesh::Real T ) const;

    libMesh::Real deriv( const libMesh::Real T ) const;

  private:

    ConstantViscosity();

    libMesh::Real _mu;

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real ConstantViscosity::operator()() const
  {
    return _mu;
  }

  inline
  libMesh::Real ConstantViscosity::operator()(AssemblyContext& /*context*/, unsigned int /*qp*/) const
  {
    return _mu;
  }

  inline
  libMesh::Real ConstantViscosity::operator()( const libMesh::Point& /*p*/,
                                               const libMesh::Real /*time*/ )
  {
    return _mu;
  }

  inline
  libMesh::Real ConstantViscosity::operator()( const libMesh::Real /*T*/ ) const
  {
    return (*this)();
  }

  inline
  libMesh::Real ConstantViscosity::deriv( const libMesh::Real /*T*/ ) const
  {
    return 0.0;
  }

} // end namespace GRINS

#endif // GRINS_CONSTANT_VISCOSITY_H
