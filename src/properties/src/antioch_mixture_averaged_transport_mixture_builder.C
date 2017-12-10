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
#include "grins/antioch_mixture_averaged_transport_mixture_builder.h"

namespace GRINS
{
  std::unique_ptr<Antioch::TransportMixture<libMesh::Real> >
  AntiochMixtureAveragedTransportMixtureBuilder::
  build_transport_mixture( const GetPot & input, const std::string & material,
                           const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    std::string transport_data_filename =
      input( "Materials/"+material+"/GasMixture/Antioch/transport_data", "default" );

    if( transport_data_filename == std::string("default") )
      transport_data_filename = Antioch::DefaultInstallFilename::transport_mixture();

    bool verbose_transport_read =
      input( "Materials/"+material+"/GasMixture/Antioch/verbose_transport_read", false );

    return std::unique_ptr<Antioch::TransportMixture<libMesh::Real> >
      ( new Antioch::TransportMixture<libMesh::Real>( chem_mix,
                                                      transport_data_filename,
                                                      verbose_transport_read,
                                                      Antioch::ParsingType::ASCII ) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
