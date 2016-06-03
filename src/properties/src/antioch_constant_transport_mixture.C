//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/antioch_constant_transport_mixture.h"

// GRINS
#include "grins/materials_parsing.h"

namespace GRINS
{
  template<typename Thermo>
  AntiochConstantTransportMixture<Thermo>::AntiochConstantTransportMixture( const GetPot& input,
                                                                            const std::string& material )
    : AntiochMixture(input,material)
  {
    libMesh::Real Le = MaterialsParsing::parse_lewis_number(input,material);
    _diffusivity.reset( new Antioch::ConstantLewisDiffusivity<libMesh::Real>(Le) );
    _mu.reset( new ConstantViscosity(input,material) );
    this->build_conductivity(input,material);

    return;
  }

  template<typename Thermo>
  AntiochConstantTransportMixture<Thermo>::~AntiochConstantTransportMixture()
  {
    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
