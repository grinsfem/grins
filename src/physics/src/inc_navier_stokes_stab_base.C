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
#include "grins/inc_navier_stokes_stab_base.h"

// GRINS
#include "grins/assembly_context.h"

namespace GRINS
{

  IncompressibleNavierStokesStabilizationBase::IncompressibleNavierStokesStabilizationBase( const std::string& physics_name, 
                                                                                            const GetPot& input )
    : IncompressibleNavierStokesBase(physics_name,input),
      _stab_helper( input )
  {
    return;
  }

  IncompressibleNavierStokesStabilizationBase::~IncompressibleNavierStokesStabilizationBase()
  {
    return;
  }

  void IncompressibleNavierStokesStabilizationBase::init_context( AssemblyContext& context )
  {
    // First call base class
    IncompressibleNavierStokesBase::init_context(context);
  
    // We need pressure derivatives
    context.get_element_fe(this->_flow_vars.p_var())->get_dphi();

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_flow_vars.u_var())->get_d2phi();

    return;
  }

  void IncompressibleNavierStokesStabilizationBase::init_variables( libMesh::FEMSystem* system )
  {
    // First call base class
    IncompressibleNavierStokesBase::init_variables(system);

    _stab_helper.init(*system);

    return;
  }

} // namespace GRINS
