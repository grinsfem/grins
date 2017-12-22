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
#include "grins/elastic_cable_abstract.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ElasticCableAbstract::ElasticCableAbstract( const PhysicsName& physics_name,
                                              const GetPot& input )
    : SolidMechanicsAbstract(physics_name,input),
      _A( 0.0 ),
      _rho(0.0)
  {
    MaterialsParsing::read_property( input,
                                     "CrossSectionalArea",
                                     PhysicsNaming::elastic_cable(),
                                     (*this),
                                     _A );

    MaterialsParsing::read_property( input,
                                     "Density",
                                     PhysicsNaming::elastic_cable(),
                                     (*this),
                                     _rho );

  }

  void ElasticCableAbstract::init_context( AssemblyContext& context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_dphidxi();

    // Need for constructing metric tensors
    this->get_fe(context)->get_dxyzdxi();
    this->get_fe(context)->get_dxidx();
    this->get_fe(context)->get_dxidy();
    this->get_fe(context)->get_dxidz();
  }

} // end namespace GRINS
