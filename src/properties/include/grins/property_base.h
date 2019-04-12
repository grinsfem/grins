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

#ifndef GRINS_PROPERTY_BASE_H
#define GRINS_PROPERTY_BASE_H

//GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/point.h"

// C++
#include <string>

namespace GRINS
{
  //! Base class for any pointwise function we might use in Physics or other places
  /*! We use a CRTP pattern for static polymorphism and to enforce an interface
      on these inhomogenous properties. */
  template<typename DerivedType>
  class PropertyBase
  {
  public:

    PropertyBase() = default;

    virtual ~PropertyBase() = default;

    libMesh::Real operator()(AssemblyContext & context, unsigned int qp) const;

    libMesh::Real operator()(const libMesh::Point & p, const libMesh::Real time);

  };

  template<typename DerivedType>
  inline
  libMesh::Real PropertyBase<DerivedType>::operator()(AssemblyContext & context, unsigned int qp) const
  {
    return static_cast<const DerivedType*>(this)->op_context_impl(context,qp);
  }

  template<typename DerivedType>
  inline
  libMesh::Real PropertyBase<DerivedType>::operator()(const libMesh::Point & p, const libMesh::Real time)
  {
    return static_cast<DerivedType*>(this)->op_point_impl(p,time);
  }

} // end namespace GRINS

#endif // GRINS_PROPERTY_BASE_H
