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


#ifndef GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H
#define GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture_builder_base.h"
#include "grins/antioch_constant_transport_mixture.h"
#include "grins/property_types.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_prandtl_conductivity.h"

// Antioch
#include "antioch/constant_lewis_diffusivity.h"

namespace GRINS
{

  class AntiochConstantTransportMixtureBuilder : public AntiochMixtureBuilderBase
  {
  public:
    AntiochConstantTransportMixtureBuilder(){}
    ~AntiochConstantTransportMixtureBuilder(){}

    template<typename KineticsThermoCurveFit,typename Conductivity>
    std::unique_ptr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
    build_mixture( const GetPot & input, const std::string & material );

    std::unique_ptr<ConstantViscosity>
    build_constant_viscosity( const GetPot & input, const std::string & material )
    {
      return std::unique_ptr<ConstantViscosity>( new ConstantViscosity(input,material) );
    }

    template<typename Conductivity>
    std::unique_ptr<Conductivity>
    build_constant_conductivity( const GetPot & input, const std::string & material )
    {
      return specialized_build_conductivity( input, material, conductivity_type<Conductivity>() );
    }

    std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> >
    build_constant_lewis_diff( const GetPot & input, const std::string & material )
    {
      libMesh::Real Le = MaterialsParsing::parse_lewis_number(input,material);
      return std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> >
        ( new Antioch::ConstantLewisDiffusivity<libMesh::Real>(Le) );
    }

  private:

    std::unique_ptr<ConstantConductivity>
    specialized_build_conductivity( const GetPot & input, const std::string & material,
                                    conductivity_type<ConstantConductivity> )
    {
      return std::unique_ptr<ConstantConductivity>( new ConstantConductivity(input,material) );
    }

    std::unique_ptr<ConstantPrandtlConductivity>
    specialized_build_conductivity( const GetPot & input, const std::string & material,
                                    conductivity_type<ConstantPrandtlConductivity> )
    {
      return std::unique_ptr<ConstantPrandtlConductivity>( new ConstantPrandtlConductivity(input,material) );
    }

  };

  template<typename KineticsThermoCurveFit,typename Conductivity>
  inline
  std::unique_ptr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
  AntiochConstantTransportMixtureBuilder::build_mixture( const GetPot & input, const std::string & material )
  {
    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix =
      this->build_chem_mix(input,material);

    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set =
      this->build_reaction_set(input,material,*chem_mix);

    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > kinetics_thermo =
      this->build_nasa_thermo_mix<KineticsThermoCurveFit>(input,material,*chem_mix);

    std::unique_ptr<ConstantViscosity> visc =
      this->build_constant_viscosity(input,material);

    std::unique_ptr<Conductivity> cond =
      this->build_constant_conductivity<Conductivity>(input,material);

    std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > diff =
      this->build_constant_lewis_diff(input,material);

    libMesh::Real min_T = this->parse_min_T(input,material);
    bool clip_negative_rho = this->parse_clip_negative_rho(input,material);

    return std::unique_ptr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
      ( new AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>
        (chem_mix, reaction_set, kinetics_thermo, visc, cond, diff, min_T, clip_negative_rho) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H
