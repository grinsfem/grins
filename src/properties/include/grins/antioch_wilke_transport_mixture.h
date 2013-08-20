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
#include "antioch/wilke_mixture.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/sutherland_parsing.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/constant_lewis_diffusivity.h"

// These are "dummy" types to help force operator overloading
namespace GRINS
{
  
}

namespace GRINS
{
  template<typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  class AntiochWilkeTransportMixture : public AntiochMixture
  {
  public:

    AntiochWilkeTransportMixture( const GetPot& input );

    virtual ~AntiochWilkeTransportMixture();

    const Antioch::WilkeMixture<libMesh::Real>& wilke_mixture() const;

    const Viscosity& viscosity() const;
    
    const Conductivity& conductivity() const;

    const Diffusivity& diffusivity() const;

    typedef AntiochChemistry ChemistryParent;
    
  protected:

    Antioch::WilkeMixture<libMesh::Real> _wilke_mixture;

    boost::scoped_ptr<Thermo> _thermo;

    boost::scoped_ptr<Viscosity> _viscosity;

    boost::scoped_ptr<Conductivity> _conductivity;

    boost::scoped_ptr<Diffusivity> _diffusivity;

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

    AntiochWilkeTransportMixture();

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

    void specialized_build_viscosity( const GetPot& /*input*/,
                                      boost::scoped_ptr<Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> > >& viscosity,
                                      viscosity_type<Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> > > )
    {
      viscosity.reset( new Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> >( *(this->_antioch_gas.get()) ) );
      
      Antioch::read_sutherland_data_ascii_default( *(viscosity.get()) );
      return;
    }

    void specialized_build_viscosity( const GetPot& /*input*/,
                                      boost::scoped_ptr<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> > >& viscosity,
                                      viscosity_type<Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> > > )
    {
      viscosity.reset( new Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> >( *(this->_antioch_gas.get()) ) );

      Antioch::read_blottner_data_ascii_default( *(viscosity.get()) );
      return;
    }

    void specialized_build_conductivity( const GetPot& /*input*/,
                                         boost::scoped_ptr<Antioch::EuckenThermalConductivity<Thermo> >& conductivity,
                                         conductivity_type<Antioch::EuckenThermalConductivity<Thermo> > )
    {
      conductivity.reset( new Antioch::EuckenThermalConductivity<Thermo>( *_thermo.get() ) );
      return;
    }

    void specialized_build_diffusivity( const GetPot& input,
                                        boost::scoped_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> >& diffusivity,
                                        diffusivity_type<Antioch::ConstantLewisDiffusivity<libMesh::Real> > )
    {
      if( !input.have_variable( "Physics/Antioch/Le" ) )
        {
          std::cerr << "Error: Must provide Lewis number for constant_lewis diffusivity model."
                    << std::endl;
           
          libmesh_error();
        }
       
      const libMesh::Real Le = input( "Physics/Antioch/Le", 0.0 );
       
      diffusivity.reset( new Antioch::ConstantLewisDiffusivity<libMesh::Real>( Le ) );
       
      return;
    }

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename T, typename V, typename C, typename D>
  inline
  const Antioch::WilkeMixture<libMesh::Real>& AntiochWilkeTransportMixture<T,V,C,D>::wilke_mixture() const
  {
    return _wilke_mixture;
  }

  template<typename T, typename Viscosity, typename C, typename D>
  inline
  const Viscosity& AntiochWilkeTransportMixture<T,Viscosity,C,D>::viscosity() const
  {
    return *_viscosity.get();
  }

  template<typename T, typename V, typename Conductivity, typename D>
  inline
  const Conductivity& AntiochWilkeTransportMixture<T,V,Conductivity,D>::conductivity() const
  {
    return *_conductivity.get();
  }

  template<typename T, typename V, typename C, typename Diffusivity>
  inline
  const Diffusivity& AntiochWilkeTransportMixture<T,V,C,Diffusivity>::diffusivity() const
  {
    return *_diffusivity.get();
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H
