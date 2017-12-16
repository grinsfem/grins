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
#include "grins/inc_navier_stokes_stab_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"

namespace GRINS
{

  template<class Mu>
  IncompressibleNavierStokesStabilizationBase<Mu>::IncompressibleNavierStokesStabilizationBase( const std::string& physics_name,
                                                                                                const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    _stab_helper( physics_name+"StabHelper", input )
  {
    return;
  }

  template<class Mu>
  IncompressibleNavierStokesStabilizationBase<Mu>::~IncompressibleNavierStokesStabilizationBase()
  {
    return;
  }

  template<class Mu>
  void IncompressibleNavierStokesStabilizationBase<Mu>::init_context( AssemblyContext& context )
  {
    // First call base class
    IncompressibleNavierStokesBase<Mu>::init_context(context);

    // We need pressure derivatives
    context.get_element_fe(this->_press_var.p())->get_dphi();

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    return;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokesStabilizationBase);
