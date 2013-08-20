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

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// This class
#include "grins/antioch_wilke_transport_mixture.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename T, typename V, typename C, typename D>
  AntiochWilkeTransportMixture<T,V,C,D>::AntiochWilkeTransportMixture( const GetPot& input )
    : AntiochMixture(input),
      _wilke_mixture(*(this->_antioch_gas.get())),
      _thermo(NULL),
      _viscosity(NULL),
      _conductivity(NULL),
      _diffusivity(NULL)
  {
    this->build_thermo( input );

    this->build_viscosity( input );

    this->build_conductivity( input );

    this->build_diffusivity( input );

    return;
  }

  template<typename T, typename V, typename C, typename D>
  AntiochWilkeTransportMixture<T,V,C,D>::~AntiochWilkeTransportMixture()
  {
    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
