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
#include "grins/elastic_cable_base.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ElasticCableBase::ElasticCableBase( const PhysicsName& physics_name,
                                      const GetPot& input )
    : Physics(physics_name,input),
      _A( 0.0 ),
      _rho(0.0),
      _disp_vars(input,physics_name)
  {
    MaterialsParsing::read_property( input,
                                     "Physics/"+physics_name+"/A",
                                     "CrossSectionalArea",
                                     elastic_cable,
                                     (*this),
                                     _A );

    MaterialsParsing::read_property( input,
                                     "Physics/"+physics_name+"/rho",
                                     "Density",
                                     elastic_cable,
                                     (*this),
                                     _rho );

  }

  ElasticCableBase::~ElasticCableBase()
  {
    return;
  }

  void ElasticCableBase::init_variables( libMesh::FEMSystem* system )
  {
    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }


  void ElasticCableBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u_var());
    system->time_evolving(_disp_vars.v_var());
    system->time_evolving(_disp_vars.w_var());

    return;
  }

  void ElasticCableBase::init_context( AssemblyContext& context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_dphidxi();

    // Need for constructing metric tensors
    this->get_fe(context)->get_dxyzdxi();
    this->get_fe(context)->get_dxidx();
    this->get_fe(context)->get_dxidy();
    this->get_fe(context)->get_dxidz();

    return;
  }

} // end namespace GRINS
