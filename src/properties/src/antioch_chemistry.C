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
#include "grins/antioch_chemistry.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/antioch_mixture_builder_base.h"

// Antioch
#include "antioch/default_filename.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  AntiochChemistry::AntiochChemistry( const GetPot& input,
                                      const std::string& material )
    : ParameterUser("AntiochChemistry")
  {
    _antioch_gas = AntiochMixtureBuilderBase().build_chem_mix(input,material);
  }

  AntiochChemistry::AntiochChemistry
  ( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture )
    : ParameterUser("AntiochChemistry")
  {
    _antioch_gas = std::move(chem_mixture);
  }

  std::string AntiochChemistry::species_name( unsigned int species_index ) const
  {
    libmesh_assert_less(species_index, _antioch_gas->n_species());

    return _antioch_gas->species_inverse_name_map().find(species_index)->second;
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
