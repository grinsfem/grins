//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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

namespace GRINS
{

  HeatTransferStabilizationBase::HeatTransferStabilizationBase( const std::string& physics_name, 
                                                                const GetPot& input )
    : HeatTransferBase(physics_name,input),
      _stab_helper( input )
  {
    this->read_input_options(input);

    return;
  }

  HeatTransferStabilizationBase::~HeatTransferStabilizationBase()
  {
    return;
  }

  void HeatTransferStabilizationBase::init_variables( libMesh::FEMSystem* system )
  {
    // First call base class
    HeatTransferBase::init_variables(system);

    _stab_helper.init(*system);

    return;
  }

  void HeatTransferStabilizationBase::init_context( AssemblyContext& context )
  {
    // First call base class
    HeatTransferBase::init_context(context);

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_temp_vars.T_var())->get_d2phi();

    return;
  }

} // namespace GRINS
