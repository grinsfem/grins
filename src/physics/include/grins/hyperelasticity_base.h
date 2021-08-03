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

#ifndef GRINS_HYPERELASTICITY_BASE_H
#define GRINS_HYPERELASTICITY_BASE_H

// GRINS
#include "grins/cartesian_solid_mechanics.h"
#include "grins/hyperelastic_strain_energy.h"

namespace GRINS
{
  template<unsigned int Dim,typename StrainEnergy>
  class HyperelasticityBase : public CartesianSolidMechanics<Dim>
  {
  public:

    HyperelasticityBase( const PhysicsName & physics_name,
                         const PhysicsName & core_physics_name,
                         const GetPot & input );

    virtual ~HyperelasticityBase() = default;

  protected:

    std::unique_ptr<HyperelasticStrainEnergy<StrainEnergy>> _strain_energy;
  };

} // end namespace GRINS

#endif // GRINS_HYPERELASTICITY_BASE_H
