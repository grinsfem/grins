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

#ifndef GRINS_SHARED_PTR_H
#define GRINS_SHARED_PTR_H

#include "grins_config.h"
#include "libmesh/libmesh_config.h"

// Here, we just piggy back on libMesh's test for shared_ptr
#ifdef LIBMESH_HAVE_CXX11_SHARED_PTR
#include <memory>
#elif GRINS_HAVE_BOOST_SHARED_PTR_HPP
#include <boost/shared_ptr.hpp>
#endif

namespace GRINS
{
#ifdef LIBMESH_HAVE_CXX11_SHARED_PTR
  template<typename T>
  class SharedPtr : public std::shared_ptr<T>
  {
  public:
    SharedPtr() : std::shared_ptr<T>() {};
    SharedPtr( T* ptr ) : std::shared_ptr<T>(ptr) {};
  };
#elif GRINS_HAVE_BOOST_SHARED_PTR_HPP
  template<typename T>
  class SharedPtr : public boost::shared_ptr<T>
  {
  public:
    SharedPtr() : boost::shared_ptr<T>() {};
    SharedPtr( T* ptr ) : boost::shared_ptr<T>(ptr) {};
  };
#else
#     error "No valid definition for shared_ptr found!"
#endif
}

#endif // GRINS_SHARED_PTR_H
