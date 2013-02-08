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

// This class
#include "grins/ideal_gas_mixture.h"

// GRINS
#include "grins/cea_thermo.h"
#include "grins/cantera_thermo.h"
#include "grins/constant_transport.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"
#include "grins/grins_kinetics.h"

namespace GRINS
{
  template<typename Thermo, typename Transport, typename Kinetics>
  IdealGasMixture<Thermo,Transport,Kinetics>::IdealGasMixture( const GetPot& input )
    : _chem_mixture( this->read_species_list(input) ), /* This *must* be done before the others */
      _thermo( input, _chem_mixture ),
      _transport( input, _chem_mixture ),
      _kinetics( input, _chem_mixture )
  {
    return;
  }

  template<typename Thermo, typename Transport, typename Kinetics>
  IdealGasMixture<Thermo,Transport,Kinetics>::~IdealGasMixture()
  {
    return;
  }

  template<typename Thermo, typename Transport, typename Kinetics>
  std::vector<std::string> IdealGasMixture<Thermo,Transport,Kinetics>::read_species_list( const GetPot& input )
  {
    unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");
    std::vector<std::string> species_names( n_species, "DIE!" );

    for( unsigned int i = 0; i < n_species; i++ )
      {
	species_names[i] = input( "Physics/Chemistry/species", "DIE!", i );
      }

    return species_names;
  }
  

  // Instantiate
  template class IdealGasMixture<CEAThermodynamics,ConstantTransport,Kinetics>;
#ifdef GRINS_HAVE_CANTERA
  template class IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics>;
  template class IdealGasMixture<CanteraThermodynamics,ConstantTransport,CanteraKinetics>;
  template class IdealGasMixture<CEAThermodynamics,ConstantTransport,CanteraKinetics>;
#endif

} // namespace GRINS
