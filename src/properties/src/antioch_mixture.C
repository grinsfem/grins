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
#include "grins/antioch_mixture.h"
#include "grins/parameter_antioch_reset.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/antioch_mixture_builder_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_multiaccessor.h"

// Antioch
#include "antioch/read_reaction_set_data.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/stat_mech_thermo.h"

namespace GRINS
{
  template <typename KineticsThermoCurveFit>
  AntiochMixture<KineticsThermoCurveFit>::AntiochMixture( const GetPot& input,
                                                          const std::string& material )
    : AntiochChemistry(input,material),
      _reaction_set( AntiochMixtureBuilderBase().build_reaction_set(input,material,(*_antioch_gas.get())) ),
      _nasa_mixture( AntiochMixtureBuilderBase().build_nasa_thermo_mix<KineticsThermoCurveFit>(input,material,(*_antioch_gas.get())) ),
      _minimum_T( AntiochMixtureBuilderBase().parse_min_T(input,material) ),
      _clip_negative_rho( AntiochMixtureBuilderBase().parse_clip_negative_rho(input,material) )
  {
    {
      std::string warning = "==============================================\n";
      warning += "WARNING: This AntiochMixture constructor is DEPREACTED!\n";
      warning += "         Prefer alternate constructor where parsing\n";
      warning += "         is done outside this class.\n";
      warning += "==============================================\n";

      libmesh_warning(warning);
    }

    this->build_stat_mech_ref_correction();
  }

  template <typename KineticsThermoCurveFit>
  AntiochMixture<KineticsThermoCurveFit>::AntiochMixture
  ( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > & nasa_mixture,
    libMesh::Real min_T,
    bool clip_negative_rho )
    : AntiochChemistry(chem_mixture),
      _minimum_T(min_T),
      _clip_negative_rho(clip_negative_rho)
  {
    _reaction_set = std::move(reaction_set);
    _nasa_mixture = std::move(nasa_mixture);

    this->build_stat_mech_ref_correction();
  }

  template <typename KineticsThermoCurveFit>
  void AntiochMixture<KineticsThermoCurveFit>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer ) const
  {
    // Use common code for any GRINS parameters
    AntiochChemistry::register_parameter(param_name, param_pointer);

    // Create a special setter/getter for any Antioch-only parameters
    if (param_name.find("Antioch") == 0) // name starts with Antioch
      param_pointer.push_back
        (ParameterAntiochReset
         (*this->_reaction_set.get(), param_name));
  }

  template <typename KineticsThermoCurveFit>
  void AntiochMixture<KineticsThermoCurveFit>::build_stat_mech_ref_correction()
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
