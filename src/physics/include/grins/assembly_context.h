//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_ASSEMBLY_CONTEXT_H
#define GRINS_ASSEMBLY_CONTEXT_H

// libMesh
#include "libmesh/fem_context.h"

// GRINS
#include "grins/cached_values.h"

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  class AssemblyContext : public libMesh::FEMContext
  {
  public:

    AssemblyContext( const libMesh::System& system );
    ~AssemblyContext();

    CachedValues & get_cached_values()
    { return _cached_values; }

    const CachedValues & get_cached_values() const
    { return _cached_values; }

    MultiphysicsSystem & get_multiphysics_system();

    const MultiphysicsSystem & get_multiphysics_system() const;

  protected:

    CachedValues _cached_values;

  };

} // end namespace GRINS

#endif // GRINS_ASSEMBLY_CONTEXT_H
