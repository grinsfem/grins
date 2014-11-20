//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/elastic_membrane_base.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ElasticMembraneBase::ElasticMembraneBase( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _disp_vars(input,physics_name)
  {
    return;
  }
  
  ElasticMembraneBase::~ElasticMembraneBase()
  {
    return;
  }
  
  void ElasticMembraneBase::init_variables( libMesh::FEMSystem* system )
  {
    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }

  void ElasticMembraneBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u_var());
    system->time_evolving(_disp_vars.v_var());
    system->time_evolving(_disp_vars.w_var());

    return;
  }

  void ElasticMembraneBase::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_disp_vars.u_var())->get_JxW();
    context.get_element_fe(_disp_vars.u_var())->get_phi();
    context.get_element_fe(_disp_vars.u_var())->get_dphidxi();
    context.get_element_fe(_disp_vars.u_var())->get_dphideta();

    // Need for constructing metric tensors
    context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();
    context.get_element_fe(_disp_vars.u_var())->get_dxyzdeta();
    context.get_element_fe(_disp_vars.u_var())->get_dxidx();
    context.get_element_fe(_disp_vars.u_var())->get_dxidy();
    context.get_element_fe(_disp_vars.u_var())->get_dxidz();
    context.get_element_fe(_disp_vars.u_var())->get_detadx();
    context.get_element_fe(_disp_vars.u_var())->get_detady();
    context.get_element_fe(_disp_vars.u_var())->get_detadz();

    return;
  }

} // end namespace GRINS
