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
#include "grins/solid_mechanics_abstract.h"

// GRINS
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  SolidMechanicsAbstract::SolidMechanicsAbstract(const PhysicsName& physics_name,
                                                 const GetPot& input )
    : Physics(physics_name,input),
      _disp_vars(input,physics_name,false,true)// is_2D = false, is_3D = true
  {}

  void SolidMechanicsAbstract::init_variables( libMesh::FEMSystem* system )
  {
    _disp_vars.init(system);
  }

  void SolidMechanicsAbstract::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u());
    system->time_evolving(_disp_vars.v());
    system->time_evolving(_disp_vars.w());
  }

} // end namespace GRINS
