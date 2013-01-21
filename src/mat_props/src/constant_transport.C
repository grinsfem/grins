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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins/constant_transport.h"

namespace GRINS
{
  ConstantTransport::ConstantTransport( const GetPot& input, const ChemicalMixture& chem_mixture )
    : _chem_mixture(chem_mixture),
      _mu( input( "Physics/Transport/viscosity", -1.0 ) ),
      _k( input( "Physics/Transport/thermal_conductivity", -1.0 ) ),
      _Le( input( "Physics/Transport/Lewis_number", -1.0 ) )
  {
    libmesh_assert_greater( _mu, 0.0 );
    libmesh_assert_greater( _k, 0.0 );
    libmesh_assert_greater( _Le, 0.0 );

    return;
  }

  ConstantTransport::~ConstantTransport()
  {
    return;
  }

} // namespace GRINS
