//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

  void ElasticMembraneBase::init_variables( libMesh::FEMSystem* system )
  {
    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }

  void ElasticMembraneBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u());
    system->time_evolving(_disp_vars.v());
    system->time_evolving(_disp_vars.w());

    return;
  }

  void ElasticMembraneBase::init_context( AssemblyContext& context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_dphidxi();
    this->get_fe(context)->get_dphideta();

    // Need for constructing metric tensors
    this->get_fe(context)->get_dxyzdxi();
    this->get_fe(context)->get_dxyzdeta();
    this->get_fe(context)->get_dxidx();
    this->get_fe(context)->get_dxidy();
    this->get_fe(context)->get_dxidz();
    this->get_fe(context)->get_detadx();
    this->get_fe(context)->get_detady();
    this->get_fe(context)->get_detadz();

    return;
  }

} // end namespace GRINS
