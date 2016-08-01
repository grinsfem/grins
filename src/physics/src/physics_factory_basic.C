//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/physics_factory_basic.h"

// GRINS
#include "grins/physics_naming.h"
#include "grins/scalar_ode.h"
#include "grins/boussinesq_buoyancy.h"
#include "grins/axisym_boussinesq_buoyancy.h"
#include "grins/elastic_membrane_constant_pressure.h"
#include "grins/constant_source_term.h"
#include "grins/parsed_source_term.h"
#include "grins/convection_diffusion.h"

namespace GRINS
{
  template<typename DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryBasic<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                  const std::string& physics_name )
  {
    return libMesh::UniquePtr<Physics>( new DerivedPhysics(physics_name,input) );
  }

  // Instantiate all the "Basic" Physics factories.
  // These shouldn't be directly used by the user, we just need to instantiate them.
  PhysicsFactoryBasic<ScalarODE> grins_factory_scalar_ode(PhysicsNaming::scalar_ode());
  PhysicsFactoryBasic<BoussinesqBuoyancy> grins_factory_boussinesq(PhysicsNaming::boussinesq_buoyancy());
  // This one needs to die. Regular Boussinesq should handle the axisymmetry
  PhysicsFactoryBasic<AxisymmetricBoussinesqBuoyancy> grins_factory_axi_boussinesq(PhysicsNaming::axisymmetric_boussinesq_buoyancy());
  PhysicsFactoryBasic<ElasticMembraneConstantPressure> grins_factory_elastic_membrane_pressure(PhysicsNaming::elastic_membrane_constant_pressure());
  // This one needs to die and just have the parsed version
  PhysicsFactoryBasic<ConstantSourceTerm> grins_factory_constant_source_term(PhysicsNaming::constant_source_term());
  PhysicsFactoryBasic<ParsedSourceTerm> grins_factory_parsed_source_term(PhysicsNaming::parsed_source_term());
  PhysicsFactoryBasic<ConvectionDiffusion> grins_factory_convection_difffusion(PhysicsNaming::convection_diffusion());

} // end namespace GRINS
