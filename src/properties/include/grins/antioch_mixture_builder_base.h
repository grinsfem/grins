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


#ifndef GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
#define GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/materials_parsing.h"
#include "grins/property_types.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_prandtl_conductivity.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/nasa_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/constant_lewis_diffusivity.h"

// libMesh
#include "libmesh/auto_ptr.h" // libMesh::UniquePtr
#include "libmesh/getpot.h"

// C++
#include <string>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Base class building Antioch mixture wrappers
  /*! This class only worries about building the kinetics
      and the thermo associated with kinetics. Subclasses
      will handle thermo and transport. */
  class AntiochMixtureBuilderBase
  {
  public:
    AntiochMixtureBuilderBase(){}
    ~AntiochMixtureBuilderBase(){}

    libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> >
    build_chem_mix( const GetPot & input, const std::string & material );

    libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> >
    build_reaction_set( const GetPot & input, const std::string & material,
                        const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    template<typename KineticsThermoCurveFit>
    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
    build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                           const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    libMesh::UniquePtr<ConstantViscosity>
    build_constant_viscosity( const GetPot & input, const std::string & material )
    {
      return libMesh::UniquePtr<ConstantViscosity>( new ConstantViscosity(input,material) );
    }

    template<typename Conductivity>
    libMesh::UniquePtr<Conductivity>
    build_constant_conductivity( const GetPot & input, const std::string & material )
    {
      return specialized_build_conductivity( input, material, conductivity_type<Conductivity>() );
    }

    libMesh::UniquePtr<Antioch::ConstantLewisDiffusivity<libMesh::Real> >
    build_constant_lewis_diff( const GetPot & input, const std::string & material )
    {
      libMesh::Real Le = MaterialsParsing::parse_lewis_number(input,material);
      return libMesh::UniquePtr<Antioch::ConstantLewisDiffusivity<libMesh::Real> >
        ( new Antioch::ConstantLewisDiffusivity<libMesh::Real>(Le) );
    }

  protected:

    void parse_nasa_data
    ( Antioch::NASAThermoMixture<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> > & nasa_mixture,
      const GetPot & input, const std::string & material)
    {
      std::string cea_data_filename = input( "Materials/"+material+"/GasMixture/Antioch/cea_data", "default" );

      if( cea_data_filename == std::string("default") )
        cea_data_filename = Antioch::DefaultInstallFilename::thermo_data();

      Antioch::read_nasa_mixture_data( nasa_mixture, cea_data_filename, Antioch::ASCII, true );
    }

    libMesh::UniquePtr<ConstantConductivity>
    specialized_build_conductivity( const GetPot & input, const std::string & material,
                                    conductivity_type<ConstantConductivity> )
    {
      return libMesh::UniquePtr<ConstantConductivity>( new ConstantConductivity(input,material) );
    }

    libMesh::UniquePtr<ConstantPrandtlConductivity>
    specialized_build_conductivity( const GetPot & input, const std::string & material,
                                    conductivity_type<ConstantPrandtlConductivity> )
    {
      return libMesh::UniquePtr<ConstantPrandtlConductivity>( new ConstantPrandtlConductivity(input,material) );
    }

  };

  template<typename KineticsThermoCurveFit>
  inline
  libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
  AntiochMixtureBuilderBase::build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                                                    const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
      nasa_mixture( new Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit>(chem_mix) );

    this->parse_nasa_data( *nasa_mixture, input, material);

    return nasa_mixture;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
