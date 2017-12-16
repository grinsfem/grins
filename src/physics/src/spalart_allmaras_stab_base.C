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
#include "grins/spalart_allmaras_stab_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/turbulence_models_macro.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmarasStabilizationBase<Mu>::SpalartAllmarasStabilizationBase( const std::string& physics_name,
                                                                          const GetPot& input )
    : SpalartAllmaras<Mu>(physics_name,input),
    _stab_helper( physics_name+"StabHelper", input )
  {}

  template<class Mu>
  void SpalartAllmarasStabilizationBase<Mu>::init_context( AssemblyContext& context )
  {
    // First call base class
    SpalartAllmaras<Mu>::init_context(context);

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_turbulence_vars.nu())->get_d2phi();
  }

  template<class Mu>
  void SpalartAllmarasStabilizationBase<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    SpalartAllmaras<Mu>::register_parameter(param_name, param_pointer);
    this->_stab_helper.register_parameter(param_name, param_pointer);
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmarasStabilizationBase);
