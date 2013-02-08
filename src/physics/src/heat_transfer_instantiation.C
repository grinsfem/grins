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
//
// $Id: heat_transfer.C 36679 2013-02-05 18:51:10Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// GRINS
#include "grins/constant_conductivity.h"

// Required source
#include "heat_transfer_base.C"
#include "heat_transfer.C"
#include "heat_transfer_stab_base.C"
#include "heat_transfer_adjoint_stab.C"

 // Instantiate
namespace GRINS
{
  template class HeatTransfer<GRINS::ConstantConductivity>;
  template class HeatTransferBase<GRINS::ConstantConductivity>;
  template class HeatTransferStabilizationBase<GRINS::ConstantConductivity>;
  template class HeatTransferAdjointStabilization<GRINS::ConstantConductivity>;
}
