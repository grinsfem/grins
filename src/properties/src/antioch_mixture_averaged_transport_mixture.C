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
#include "grins/antioch_mixture_averaged_transport_mixture.h"

// GRINS
#include "grins/antioch_mixture_averaged_transport_mixture_builder.h"

// Antioch
#include "antioch/default_filename.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename KT, typename T, typename V, typename C, typename D>
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::AntiochMixtureAveragedTransportMixture
  ( const GetPot& input,const std::string& material )
    : AntiochMixture<KT>(input,material)
  {
    {
      std::string warning = "==============================================\n";
      warning += "WARNING: This AntiochMixtureAveragedTransportMixture constructor is DEPREACTED!\n";
      warning += "         Prefer alternate constructor where parsing\n";
      warning += "         is done outside this class.\n";
      warning += "==============================================\n";

      libmesh_warning(warning);
    }

    AntiochMixtureAveragedTransportMixtureBuilder builder;

    _trans_mixture = builder.build_transport_mixture(  input, material, *(this->_antioch_gas.get()) );

    _wilke_mixture = builder.build_mix_avg_trans_mixture(*(this->_trans_mixture.get()));

    _thermo = builder.build_gas_thermo<KT,T>( *(this->_antioch_gas.get()), *(this->_nasa_mixture.get()) );

    _viscosity = builder.build_viscosity<V>(input, material, (*(this->_trans_mixture.get())));

    _conductivity = builder.build_conductivity<C>(*(this->_trans_mixture.get()), *(this->_thermo.get()));

    _diffusivity = builder.build_diffusivity<D>(input, material, (*(this->_trans_mixture.get())));
  }

  template<typename KT, typename T, typename V, typename C, typename D>
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::AntiochMixtureAveragedTransportMixture
  ( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KT> > & kinetics_thermo_mix,
    std::unique_ptr<T> & gas_thermo,
    std::unique_ptr<Antioch::TransportMixture<libMesh::Real> > & trans_mix,
    std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> > & wilke_mix,
    std::unique_ptr<Antioch::MixtureViscosity<V,libMesh::Real> > & visc,
    std::unique_ptr<Antioch::MixtureConductivity<C,libMesh::Real> > & conductivity,
    std::unique_ptr<Antioch::MixtureDiffusion<D,libMesh::Real> > & diff,
    libMesh::Real min_T,
    bool clip_negative_rho )
  : AntiochMixture<KT>(chem_mixture,reaction_set,kinetics_thermo_mix,min_T,clip_negative_rho),
    _trans_mixture(std::move(trans_mix)),
    _wilke_mixture(std::move(wilke_mix)),
    _thermo(std::move(gas_thermo)),
    _viscosity(std::move(visc)),
    _conductivity(std::move(conductivity)),
    _diffusivity(std::move(diff))
  {}

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
