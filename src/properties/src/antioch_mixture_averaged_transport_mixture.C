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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// This class
#include "grins/antioch_mixture_averaged_transport_mixture.h"

// Antioch
#include "antioch/default_filename.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename T, typename V, typename C, typename D>
  AntiochMixtureAveragedTransportMixture<T,V,C,D>::AntiochMixtureAveragedTransportMixture( const GetPot& input )
    : AntiochMixture(input),
      _trans_mixture(NULL),
      _wilke_mixture(NULL),
      _thermo(NULL),
      _viscosity(NULL),
      _conductivity(NULL),
      _diffusivity(NULL)
  {
    std::string transport_data_filename = input( "Physics/Antioch/transport_data", "default" );
    if( transport_data_filename == std::string("default") )
      transport_data_filename = Antioch::DefaultInstallFilename::transport_mixture();

    bool verbose_transport_read = input( "Physics/Antioch/verbose_transport_read", false );

    _trans_mixture.reset( new Antioch::TransportMixture<libMesh::Real>( *(_antioch_gas.get()),
                                                                        transport_data_filename,
                                                                        verbose_transport_read,
                                                                        Antioch::ParsingType::ASCII ) );

    _wilke_mixture.reset( new Antioch::MixtureAveragedTransportMixture<libMesh::Real>(*(_trans_mixture.get()) ) );

    this->build_thermo( input );

    this->build_viscosity( input );

    this->build_conductivity( input );

    this->build_diffusivity( input );

    return;
  }

  template<typename T, typename V, typename C, typename D>
  AntiochMixtureAveragedTransportMixture<T,V,C,D>::~AntiochMixtureAveragedTransportMixture()
  {
    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
