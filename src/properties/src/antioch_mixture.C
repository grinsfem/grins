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
#include "grins/antioch_mixture.h"
#include "grins/parameter_antioch_reset.h"

// GRINS
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_multiaccessor.h"

// Antioch
#include "antioch/read_reaction_set_data.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/stat_mech_thermo.h"

namespace GRINS
{
  AntiochMixture::AntiochMixture( const GetPot& input,
                                  const std::string& material )
    : AntiochChemistry(input,material),
      _reaction_set( new Antioch::ReactionSet<libMesh::Real>( (*_antioch_gas.get()) ) ),
      _cea_mixture( new Antioch::CEAThermoMixture<libMesh::Real>( (*_antioch_gas.get()) ) ),
      _minimum_T( input( "Materials/"+material+"/GasMixture/Antioch/minimum_T", -std::numeric_limits<libMesh::Real>::max() ) ),
      _clip_negative_rho( input( "Materials/"+material+"/GasMixture/Antioch/clip_negative_rho", false) )
  {
    std::string kinetics_data_filename = MaterialsParsing::parse_chemical_kinetics_datafile_name( input, material );

    bool verbose_read = input("screen-options/verbose_kinetics_read", false );

    Antioch::read_reaction_set_data_xml<libMesh::Real>( kinetics_data_filename, verbose_read, *_reaction_set.get() );

    std::string cea_data_filename = input( "Materials/"+material+"/GasMixture/Antioch/cea_data", "default" );
    if( cea_data_filename == std::string("default") )
      cea_data_filename = Antioch::DefaultInstallFilename::thermo_data();

    Antioch::read_cea_mixture_data_ascii( *_cea_mixture.get(), cea_data_filename );

    this->build_stat_mech_ref_correction();
  }

  void AntiochMixture::register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    // Use common code for any GRINS parameters
    AntiochChemistry::register_parameter(param_name, param_pointer);

    // Create a special setter/getter for any Antioch-only parameters
    if (param_name.find("Antioch") == 0) // name starts with Antioch
      param_pointer.push_back
        (ParameterAntiochReset
          (*this->_reaction_set.get(), param_name));
  }


  void AntiochMixture::build_stat_mech_ref_correction()
  {
    Antioch::StatMechThermodynamics<libMesh::Real> thermo( *(this->_antioch_gas.get()) );

    _h_stat_mech_ref_correction.resize(this->n_species());

    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        _h_stat_mech_ref_correction[s] = -thermo.h_tot( s, 298.15 ) + thermo.e_0(s);
      }
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
