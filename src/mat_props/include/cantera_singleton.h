//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_CANTERA_SINGLETON_H
#define GRINS_CANTERA_SINGLETON_H

// GRINS
#include "grins_config.h"

// Boost
#include "boost/tr1/memory.hpp"

// Cantera
#ifdef GRINS_HAVE_CANTERA
#include "cantera/IdealGasMix.h"
#endif

// libMesh
#include "libmesh/getpot.h"

namespace
{
  static libMesh::Threads::spin_mutex cantera_mutex;
}

namespace GRINS
{
#ifdef GRINS_HAVE_CANTERA

  class CanteraSingleton
  {
  public:

    static Cantera::IdealGasMix& cantera_instance( const GetPot& input );

    static Cantera::IdealGasMix& cantera_instance();

  private:

    CanteraSingleton();
    
    //! Pointer to Cantera IdealGasMix instance.
    /*! Use a smart pointer here so that it gets properly destroyed at the end
      of the program.
      \todo Use a lighter-weight smart pointer. Thinking std::unique_ptr when it's ubiquitous enough */
    static std::tr1::shared_ptr<Cantera::IdealGasMix> _cantera;
  };

#endif //GRINS_HAVE_CANTERA

} // namespace GRINS

#endif //GRINS_CANTERA_SINGLETON_H
