//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H
#define GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/property_types.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/getpot.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/cea_evaluator.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/mixture_averaged_transport_mixture.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/mixture_conductivity.h"
#include "antioch/mixture_diffusion.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/constant_lewis_diffusivity.h"

#include "antioch/sutherland_parsing.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/constant_lewis_diffusivity_building.h"

#ifdef ANTIOCH_HAVE_GSL

#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/gsl_spliner.h"
#include "antioch/kinetics_theory_viscosity_building.h"
#include "antioch/kinetics_theory_thermal_conductivity_building.h"

#endif // ANTIOCH_HAVE_GSL

namespace GRINS
{
  //! Wrapper class for storing state for computing Wilke transport properties using Antioch
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  class AntiochMixtureAveragedTransportMixture : public AntiochMixture
  {
  public:

    AntiochMixtureAveragedTransportMixture( const GetPot& input );

    virtual ~AntiochMixtureAveragedTransportMixture();

    const Antioch::MixtureAveragedTransportMixture<libMesh::Real>& wilke_mixture() const;

    const Antioch::MixtureViscosity<Viscosity,libMesh::Real>& viscosity() const;

    const Antioch::MixtureConductivity<Conductivity,libMesh::Real>& conductivity() const;

    const Antioch::MixtureDiffusion<Diffusivity,libMesh::Real>& diffusivity() const;

    typedef AntiochChemistry ChemistryParent;

  protected:

    boost::scoped_ptr<Antioch::TransportMixture<libMesh::Real> > _trans_mixture;

    boost::scoped_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> > _wilke_mixture;

    boost::scoped_ptr<Thermo> _thermo;

    boost::scoped_ptr<Antioch::MixtureViscosity<Viscosity,libMesh::Real> > _viscosity;

    boost::scoped_ptr<Antioch::MixtureConductivity<Conductivity,libMesh::Real> > _conductivity;

    boost::scoped_ptr<Antioch::MixtureDiffusion<Diffusivity,libMesh::Real> > _diffusivity;

    /* Below we will specialize the specialized_build_* functions to the appropriate type.
       This way, we can control how the cached transport objects get constructed
       based on the template type. This is achieved by the dummy types forcing operator
       overloading for each of the specialized types. */
    void build_thermo( const GetPot& input )
    { specialized_build_thermo( input, _thermo, thermo_type<Thermo>() ); }

    void build_viscosity( const GetPot& input )
    { specialized_build_viscosity( input, _viscosity, viscosity_type<Viscosity>() ); }

    void build_conductivity( const GetPot& input )
    { specialized_build_conductivity( input, _conductivity, conductivity_type<Conductivity>() ); }

    void build_diffusivity( const GetPot& input )
    { specialized_build_diffusivity( input, _diffusivity, diffusivity_type<Diffusivity>() ); }

  private:

    AntiochMixtureAveragedTransportMixture();

    void specialized_build_thermo( const GetPot& /*input*/,
                                   boost::scoped_ptr<Antioch::StatMechThermodynamics<libMesh::Real> >& thermo,
                                   thermo_type<Antioch::StatMechThermodynamics<libMesh::Real> > )
    {
      thermo.reset( new Antioch::StatMechThermodynamics<libMesh::Real>( *(this->_antioch_gas.get()) ) );
      return;
    }

    void specialized_build_thermo( const GetPot& /*input*/,
                                   boost::scoped_ptr<Antioch::CEAEvaluator<libMesh::Real> >& thermo,
                                   thermo_type<Antioch::CEAEvaluator<libMesh::Real> > )
    {
      thermo.reset( new Antioch::CEAEvaluator<libMesh::Real>( this->cea_mixture() ) );
      return;
    }

    void specialized_build_viscosity( const GetPot& input,
                                      boost::scoped_ptr<Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real>,libMesh::Real> >& viscosity,
                                      viscosity_type<Antioch::SutherlandViscosity<libMesh::Real> > )
    {
      viscosity.reset( new Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real>,libMesh::Real>(*(_trans_mixture.get())) );

      std::string sutherland_data = input("Physics/Antioch/sutherland_data", "default");
      if( sutherland_data == "default" )
        sutherland_data = Antioch::DefaultInstallFilename::sutherland_data();

      Antioch::read_sutherland_data_ascii( *(viscosity.get()), sutherland_data );
    }

    void specialized_build_viscosity( const GetPot& input,
                                      boost::scoped_ptr<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real>,libMesh::Real> >& viscosity,
                                      viscosity_type<Antioch::BlottnerViscosity<libMesh::Real> > )
    {
      viscosity.reset( new Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real>,libMesh::Real>(*(_trans_mixture.get())) );

      std::string blottner_data = input("Physics/Antioch/blottner_data", "default");
      if( blottner_data == "default" )
        blottner_data = Antioch::DefaultInstallFilename::blottner_data();

      Antioch::read_blottner_data_ascii( *(viscosity.get()), blottner_data );
    }

#ifdef ANTIOCH_HAVE_GSL
    void specialized_build_viscosity( const GetPot& /*input*/,
                                      boost::scoped_ptr<Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >& viscosity,
                                      viscosity_type<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner> > )
    {
      viscosity.reset( new Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real>(*(_trans_mixture.get())) );

      Antioch::build_kinetics_theory_viscosity<libMesh::Real,Antioch::GSLSpliner>( *(viscosity.get()) );
    }
#endif // ANTIOCH_HAVE_GSL

    void specialized_build_conductivity( const GetPot& /*input*/,
                                         boost::scoped_ptr<Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<Thermo>,libMesh::Real> >& conductivity,
                                         conductivity_type<Antioch::EuckenThermalConductivity<Thermo> > )
    {
      conductivity.reset( new Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<Thermo>,libMesh::Real>(*(_trans_mixture.get())) );
      Antioch::build_eucken_thermal_conductivity<Thermo,libMesh::Real>(*(conductivity.get()),*(_thermo.get()));
      return;
    }

#ifdef ANTIOCH_HAVE_GSL
    void specialized_build_conductivity( const GetPot& /*input*/,
                                         boost::scoped_ptr<Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,libMesh::Real> >& conductivity,
                                         conductivity_type<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real> > )
    {
      conductivity.reset( new Antioch::MixtureConductivity<Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,libMesh::Real>(*(_trans_mixture.get())) );

      Antioch::build_kinetics_theory_thermal_conductivity<Thermo,libMesh::Real>( *(conductivity.get()), *(_thermo.get()) );
    }
#endif // ANTIOCH_HAVE_GSL

    void specialized_build_diffusivity( const GetPot& input,
                                        boost::scoped_ptr<Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<libMesh::Real>,libMesh::Real> >& diffusivity,
                                        diffusivity_type<Antioch::ConstantLewisDiffusivity<libMesh::Real> > )
    {
      if( !input.have_variable( "Physics/Antioch/Le" ) )
        {
          std::cerr << "Error: Must provide Lewis number for constant_lewis diffusivity model."
                    << std::endl;

          libmesh_error();
        }

      const libMesh::Real Le = input( "Physics/Antioch/Le", 0.0 );

      diffusivity.reset( new Antioch::MixtureDiffusion<Antioch::ConstantLewisDiffusivity<libMesh::Real>,libMesh::Real>(*(_trans_mixture.get())) );

      Antioch::build_constant_lewis_diffusivity<libMesh::Real>( *(diffusivity.get()), Le);
      return;
    }

#ifdef ANTIOCH_HAVE_GSL
    void specialized_build_diffusivity( const GetPot& /*input*/,
                                        boost::scoped_ptr<Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real> >& diffusivity,
                                        diffusivity_type<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > )
    {
      diffusivity.reset( new Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner>,libMesh::Real>(*(_trans_mixture.get())) );
    }
#endif // ANTIOCH_HAVE_GSL

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename T, typename V, typename C, typename D>
  inline
  const Antioch::MixtureAveragedTransportMixture<libMesh::Real>& AntiochMixtureAveragedTransportMixture<T,V,C,D>::wilke_mixture() const
  {
    return *(_wilke_mixture.get());
  }

  template<typename T, typename Viscosity, typename C, typename D>
  inline
  const Antioch::MixtureViscosity<Viscosity,libMesh::Real>& AntiochMixtureAveragedTransportMixture<T,Viscosity,C,D>::viscosity() const
  {
    return *_viscosity.get();
  }

  template<typename T, typename V, typename Conductivity, typename D>
  inline
  const Antioch::MixtureConductivity<Conductivity,libMesh::Real>& AntiochMixtureAveragedTransportMixture<T,V,Conductivity,D>::conductivity() const
  {
    return *_conductivity.get();
  }

  template<typename T, typename V, typename C, typename Diffusivity>
  inline
  const Antioch::MixtureDiffusion<Diffusivity,libMesh::Real>& AntiochMixtureAveragedTransportMixture<T,V,C,Diffusivity>::diffusivity() const
  {
    return *_diffusivity.get();
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H
