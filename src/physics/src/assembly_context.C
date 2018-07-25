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

// This class
#include "grins/assembly_context.h"

// GRINS
#include "grins/multiphysics_sys.h"

namespace GRINS
{
  AssemblyContext::AssemblyContext( const libMesh::System& system )
    : libMesh::FEMContext(system)
  {
    return;
  }

  AssemblyContext::~AssemblyContext()
  {
    return;
  }

  MultiphysicsSystem & AssemblyContext::get_multiphysics_system()
  {
    libMesh::System & base_system = const_cast<libMesh::System &>(this->get_system());

    MultiphysicsSystem & multiphysics_system =
      libMesh::cast_ref<MultiphysicsSystem &>( base_system );

    return multiphysics_system;
  }

  const MultiphysicsSystem & AssemblyContext::get_multiphysics_system() const
  {
    const MultiphysicsSystem & multiphysics_system =
      libMesh::cast_ref<const MultiphysicsSystem &>( this->get_system() );

    return multiphysics_system;
  }

} // end namespace GRINS
