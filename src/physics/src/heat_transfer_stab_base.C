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
#include "grins/heat_transfer_stab_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_conductivity.h"
#include "grins/parsed_conductivity.h"
#include "grins/heat_transfer_macros.h"

namespace GRINS
{
  template<class K>
  HeatTransferStabilizationBase<K>::HeatTransferStabilizationBase( const std::string& physics_name,
                                                                   const GetPot& input )
    : HeatTransferBase<K>(physics_name,PhysicsNaming::heat_transfer(),input),
    _stab_helper(physics_name+"StabHelper", input)
  {}

  template<class K>
  void HeatTransferStabilizationBase<K>::init_context( AssemblyContext& context )
  {
    // First call base class
    HeatTransferBase<K>::init_context(context);

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_temp_vars.T())->get_d2phi();

    return;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_HEAT_TRANSFER_SUBCLASS(HeatTransferStabilizationBase);
