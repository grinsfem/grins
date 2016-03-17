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

// This class
#include "grins/physics_factory_base.h"

// Full specialization for the Factory<Physics>
namespace libMesh
{
  template<>
  std::map<std::string, libMesh::Factory<GRINS::Physics>*>&
  libMesh::Factory<GRINS::Physics>::factory_map()
  {
    static std::map<std::string, libMesh::Factory<GRINS::Physics>*> _map;
    return _map;
  }
} // end namespace libMesh

// Definition of static members
namespace GRINS
{
  std::string PhysicsFactoryBase::_physics_name = std::string("DIE!");
  const GetPot* PhysicsFactoryBase::_input = NULL;
} // end namespace GRINS
