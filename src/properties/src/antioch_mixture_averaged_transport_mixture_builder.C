//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
    std::unique_ptr<Antioch::TransportMixture<libMesh::Real> > trans_mix;

    Antioch::ParsingType parsing_type = this->get_antioch_parsing_type(input,material);

    switch(parsing_type)
      {
      case(Antioch::ASCII):
          {
            std::string prefix(this->antioch_prefix(material));

            bool verbose_transport_read =
              input( prefix+"/verbose_transport_read", false );

            std::string transport_data_filename =
              input( prefix+"/transport_data", "default" );

            if( transport_data_filename == std::string("default") )
              transport_data_filename = Antioch::DefaultInstallFilename::transport_mixture();

            trans_mix.reset( new Antioch::TransportMixture<libMesh::Real>( chem_mix,
                                                                           transport_data_filename,
                                                                           verbose_transport_read,
                                                                           Antioch::ParsingType::ASCII ) );

            break;
          }
      case(Antioch::XML):
        {
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          trans_mix.reset( new Antioch::TransportMixture<libMesh::Real>(chem_mix, parser.get()) );
          break;
        }
      case(Antioch::CHEMKIN):
        {
          libmesh_not_implemented();
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid Antioch parsing type!");
      }

    return trans_mix;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
