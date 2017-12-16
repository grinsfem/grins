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
#include "grins/turbulence_models_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/physics_naming.h"
#include "grins/turbulence_models_macro.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<class Mu>
  TurbulenceModelsBase<Mu>::TurbulenceModelsBase(const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input),
      _rho(input("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho", 1.0)),
      _mu(input,MaterialsParsing::material_name(input,PhysicsNaming::incompressible_navier_stokes()))
  {
    this->set_parameter(this->_rho, input, "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho", this->_rho);
  }

  template<class Mu>
  void TurbulenceModelsBase<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _mu.register_parameter(param_name, param_pointer);
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(TurbulenceModelsBase);
