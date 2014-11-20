//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/constant_conductivity.h"

//GRINS
#include "grins/grins_physics_names.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ConstantConductivity::ConstantConductivity( const GetPot& input )
    : _k( input("Materials/Conductivity/k", 0.0) )
  {
    if( !input.have_variable("Materials/Conductivity/k") )
      {
        libmesh_warning("No Materials/Conductivity/k specified!\n");

	// Try and get the conductivity from other specifications
	_k = input("Physics/"+incompressible_navier_stokes+"/k", 1.0);
      }
    return;
  }

  ConstantConductivity::~ConstantConductivity()
  {
    return;
  }

} // end namespace GRINS
