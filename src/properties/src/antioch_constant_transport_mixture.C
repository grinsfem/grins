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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

#include "grins/antioch_constant_transport_mixture.h"

namespace GRINS
{
  template<typename Thermo>
  AntiochConstantTransportMixture<Thermo>::AntiochConstantTransportMixture( const GetPot& input )
    : AntiochMixture(input),
      _mu( input("Materials/Viscosity/mu", 0.0) ),
      _conductivity(NULL),
      _diffusivity( new Antioch::ConstantLewisDiffusivity<libMesh::Real>( input("Physics/Antioch/Le", 0.0) ) )
  {
    if( !input.have_variable("Materials/Viscosity/mu") )
      {
        std::cerr << "Error: Must specify viscosity value for constant viscosity model!" << std::endl;
        libmesh_error();
      }

    if( !input.have_variable("Physics/Antioch/Le") )
      {
        std::cerr << "Error: Must provide Lewis number for constant_lewis diffusivity model."
                  << std::endl;
        libmesh_error();
      }

    this->build_conductivity(input);

    return;
  }

  template<typename Thermo>
  AntiochConstantTransportMixture<Thermo>::~AntiochConstantTransportMixture()
  {
    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
