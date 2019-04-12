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

#ifndef GRINS_ANTIOCH_MIXTURE_AVERAGED_TRANSPORT_MIXTURE_BUILDER_H
#define GRINS_ANTIOCH_MIXTURE_AVERAGED_TRANSPORT_MIXTURE_BUILDER_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture_builder_base.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"

namespace GRINS
{

  class AntiochMixtureAveragedTransportMixtureBuilder : public AntiochMixtureBuilderBase
  {
  public:
    AntiochMixtureAveragedTransportMixtureBuilder(){}
    ~AntiochMixtureAveragedTransportMixtureBuilder(){}

    template<typename KineticsThermoCurveFit, typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
    std::unique_ptr<AntiochMixtureAveragedTransportMixture<KineticsThermoCurveFit,Thermo,Viscosity,Conductivity,Diffusivity> >
    build_mixture( const GetPot & input, const std::string & material );

    std::unique_ptr<Antioch::TransportMixture<libMesh::Real> >
    build_transport_mixture( const GetPot& input, const std::string& material,
                             const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> >
    build_mix_avg_trans_mixture( const Antioch::TransportMixture<libMesh::Real> & trans_mix )
    { return std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> >
        ( new Antioch::MixtureAveragedTransportMixture<libMesh::Real>(trans_mix) ); }

    template<typename Viscosity>
    std::unique_ptr<Antioch::MixtureViscosity<Viscosity,libMesh::Real> >
    build_viscosity( const GetPot& input, const std::string& material,
                     const Antioch::TransportMixture<libMesh::Real> & trans_mix )
    { return specialized_build_viscosity( input, material, trans_mix, viscosity_type<Viscosity>() ); }

    template<typename Diffusivity>
    std::unique_ptr<Antioch::MixtureDiffusion<Diffusivity,libMesh::Real> >
    build_diffusivity( const GetPot& input, const std::string& material,
                       const Antioch::TransportMixture<libMesh::Real> & trans_mix )
    { return specialized_build_diffusivity( input, material, trans_mix, diffusivity_type<Diffusivity>() ); }

    template<typename Conductivity, typename Thermo>
    std::unique_ptr<Antioch::MixtureConductivity<Conductivity,libMesh::Real> >
    build_conductivity( const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                        const Thermo & thermo )
    { return specialized_build_conductivity<Thermo>(trans_mix, thermo, conductivity_type<Conductivity>() ); }

  private:

    std::unique_ptr<Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real>,libMesh::Real> >
    specialized_build_viscosity( const GetPot& input,
                                 const std::string& material,
                                 const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                 viscosity_type<Antioch::SutherlandViscosity<libMesh::Real> > )
    {
      std::unique_ptr<Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real>,libMesh::Real> >
        viscosity( new Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real>,libMesh::Real>(trans_mix) );

      std::string sutherland_data = input("Materials/"+material+"/GasMixture/Antioch/sutherland_data", "default");
      if( sutherland_data == "default" )
        sutherland_data = Antioch::DefaultInstallFilename::sutherland_data();

      Antioch::read_sutherland_data_ascii( *(viscosity.get()), sutherland_data );

      return viscosity;
    }

    std::unique_ptr<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real>,libMesh::Real> >
    specialized_build_viscosity( const GetPot& input,
                                 const std::string& material,
                                 const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                 viscosity_type<Antioch::BlottnerViscosity<libMesh::Real> > )
    {
      std::unique_ptr<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real>,libMesh::Real> >
        viscosity( new Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real>,libMesh::Real>(trans_mix) );

      std::string blottner_data = input("Materials/"+material+"/GasMixture/Antioch/blottner_data", "default");
      if( blottner_data == "default" )
        blottner_data = Antioch::DefaultInstallFilename::blottner_data();

      Antioch::read_blottner_data_ascii( *(viscosity.get()), blottner_data );

      return viscosity;
    }

#ifdef ANTIOCH_HAVE_GSL
    std::unique_ptr<Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >
    specialized_build_viscosity( const GetPot& /*input*/,
                                 const std::string& /*material*/,
                                 const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                 viscosity_type<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner> > )
    {
      std::unique_ptr<Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >
        viscosity( new Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real>(trans_mix) );

      Antioch::build_kinetics_theory_viscosity<libMesh::Real,Antioch::GSLSpliner>( *(viscosity.get()) );

      return viscosity;
    }
#endif // ANTIOCH_HAVE_GSL


    std::unique_ptr<Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<libMesh::Real>,libMesh::Real> >
    specialized_build_diffusivity( const GetPot& input,
                                   const std::string& material,
                                   const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                   diffusivity_type<Antioch::ConstantLewisDiffusivity<libMesh::Real> > )
    {
      libMesh::Real Le = MaterialsParsing::parse_lewis_number(input,material);

      std::unique_ptr<Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<libMesh::Real>,libMesh::Real> >
        diffusivity( new Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<libMesh::Real>,libMesh::Real>(trans_mix) );

      Antioch::build_constant_lewis_diffusivity<libMesh::Real>( *(diffusivity.get()), Le);

      return diffusivity;
    }

#ifdef ANTIOCH_HAVE_GSL
    std::unique_ptr<Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >
    specialized_build_diffusivity( const GetPot& /*input*/,
                                   const std::string& /*material*/,
                                   const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                   diffusivity_type<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > )
    {
      return std::unique_ptr<Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >
        ( new Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real>(trans_mix) );
    }
#endif // ANTIOCH_HAVE_GSL

    template<typename Thermo>
    std::unique_ptr<Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<Thermo>,libMesh::Real> >
    specialized_build_conductivity(const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                   const Thermo & thermo,
                                   conductivity_type<Antioch::EuckenThermalConductivity<Thermo> > )
    {
      std::unique_ptr<Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<Thermo>,libMesh::Real> >
        conductivity( new Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<Thermo>,libMesh::Real>(trans_mix) );
      Antioch::build_eucken_thermal_conductivity<Thermo,libMesh::Real>(*(conductivity.get()),thermo);

      return conductivity;
    }

#ifdef ANTIOCH_HAVE_GSL
    template<typename Thermo>
    std::unique_ptr<Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,libMesh::Real> >
    specialized_build_conductivity(const Antioch::TransportMixture<libMesh::Real> & trans_mix,
                                   const Thermo & thermo,
                                   conductivity_type<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real> > )
    {
      std::unique_ptr<Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,libMesh::Real> >
      conductivity( new Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,libMesh::Real>(trans_mix));

      Antioch::build_kinetics_theory_thermal_conductivity<Thermo,libMesh::Real>(*(conductivity.get()), thermo);

      return conductivity;
    }
#endif // ANTIOCH_HAVE_GSL

  };


  template<typename KT, typename T, typename V, typename C, typename D>
  inline
  std::unique_ptr<AntiochMixtureAveragedTransportMixture<KT,T,V,C,D> >
  AntiochMixtureAveragedTransportMixtureBuilder::
  build_mixture( const GetPot & input, const std::string & material )
  {
    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix =
      this->build_chem_mix(input,material);

    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set =
      this->build_reaction_set(input,material,*chem_mix);

    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KT> > nasa_mix =
      this->build_nasa_thermo_mix<KT>(input,material,*chem_mix);

    std::unique_ptr<T> gas_thermo =
      this->build_gas_thermo<KT,T>( *chem_mix, *nasa_mix );

    std::unique_ptr<Antioch::TransportMixture<libMesh::Real> > trans_mix =
      this->build_transport_mixture(input,material,*chem_mix);

    std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> > wilke_mix =
      this->build_mix_avg_trans_mixture(*trans_mix);

    std::unique_ptr<Antioch::MixtureViscosity<V,libMesh::Real> > visc =
      this->build_viscosity<V>(input,material,*trans_mix);

    std::unique_ptr<Antioch::MixtureConductivity<C,libMesh::Real> > cond =
      this->build_conductivity<C>(*trans_mix, *gas_thermo);

    std::unique_ptr<Antioch::MixtureDiffusion<D,libMesh::Real> > diff =
      this->build_diffusivity<D>(input,material,*trans_mix);

    libMesh::Real min_T = this->parse_min_T(input,material);
    bool clip_negative_rho = this->parse_clip_negative_rho(input,material);

    return std::unique_ptr<AntiochMixtureAveragedTransportMixture<KT,T,V,C,D> >
      ( new AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>
        (chem_mix, reaction_set, nasa_mix, gas_thermo,
         trans_mix, wilke_mix, visc, cond, diff,
         min_T, clip_negative_rho) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_AVERAGED_TRANSPORT_MIXTURE_BUILDER_H
