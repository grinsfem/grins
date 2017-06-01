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


#ifndef GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H
#define GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture_builder_base.h"
#include "grins/antioch_constant_transport_mixture.h"

namespace GRINS
{

  class AntiochConstantTransportMixtureBuilder : public AntiochMixtureBuilderBase
  {
  public:
    AntiochConstantTransportMixtureBuilder(){}
    ~AntiochConstantTransportMixtureBuilder(){}

    template<typename KineticsThermoCurveFit,typename Conductivity>
    libMesh::UniquePtr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
    build_mixture( const GetPot & input, const std::string & material );

  };

  template<typename KineticsThermoCurveFit,typename Conductivity>
  inline
  libMesh::UniquePtr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
  AntiochConstantTransportMixtureBuilder::build_mixture( const GetPot & input, const std::string & material )
  {
    libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix =
      this->build_chem_mix(input,material);

    libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> > reaction_set =
      this->build_reaction_set(input,material,*chem_mix);

    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > kinetics_thermo =
      this->build_nasa_thermo_mix<KineticsThermoCurveFit>(input,material,*chem_mix);

    libMesh::UniquePtr<ConstantViscosity> visc =
      this->build_constant_viscosity(input,material);

    libMesh::UniquePtr<Conductivity> cond =
      this->build_constant_conductivity<Conductivity>(input,material);

    libMesh::UniquePtr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > diff =
      this->build_constant_lewis_diff(input,material);

    libMesh::Real min_T = this->parse_min_T(input,material);
    bool clip_negative_rho = this->parse_clip_negative_rho(input,material);

    return libMesh::UniquePtr<AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity> >
      ( new AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>
        (chem_mix, reaction_set, kinetics_thermo, visc, cond, diff, min_T, clip_negative_rho) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_BUILDER_H
