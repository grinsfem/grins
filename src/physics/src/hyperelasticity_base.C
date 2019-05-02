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

// This class
#include "grins/hyperelasticity_base.h"

// GRINS
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<unsigned int Dim,typename StrainEnergy>
  HyperelasticityBase<Dim,StrainEnergy>::HyperelasticityBase
  ( const PhysicsName & physics_name, const PhysicsName & core_physics_name, const GetPot & input )
    : CartesianSolidMechanics<Dim>(physics_name,core_physics_name,input),
      _strain_energy(nullptr)
  {
    const std::string material =
      MaterialsParsing::material_name(input,core_physics_name);

    _strain_energy.reset(new StrainEnergy(input,material));
  }

} // end namespace GRINS
