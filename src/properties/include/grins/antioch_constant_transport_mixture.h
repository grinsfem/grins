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


#ifndef GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_H
#define GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// Antioch
#include "antioch/constant_lewis_diffusivity.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/property_types.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_prandtl_conductivity.h"
#include "grins/constant_viscosity.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  //! Wrapper class for storing state for constant transport properties, including Antioch::ConstantLewisDiffusivity
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit,typename Conductivity>
  class AntiochConstantTransportMixture : public AntiochMixture<KineticsThermoCurveFit>
  {
  public:

    //! Deprecated Constructor
    AntiochConstantTransportMixture( const GetPot & input, const std::string & material );

    //! Constructor with user-built objects
    AntiochConstantTransportMixture( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
                                     std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
                                     std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > & nasa_mixture,
                                     std::unique_ptr<ConstantViscosity> & visc,
                                     std::unique_ptr<Conductivity> & cond,
                                     std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > & diff,
                                     libMesh::Real min_T = -std::numeric_limits<libMesh::Real>::max(),
                                     bool clip_negative_rho = false );

    virtual ~AntiochConstantTransportMixture(){}

    libMesh::Real mu() const;

    const Conductivity& conductivity() const;

    const Antioch::ConstantLewisDiffusivity<libMesh::Real>& diffusivity() const;

    typedef AntiochChemistry ChemistryParent;

  protected:

    //! Viscosity
    /*! \todo Template on viscosity model, as we did for conductivity,
      to support parsed versions. This is going to require we update the
      API for the transport wrappers. */
    std::unique_ptr<ConstantViscosity> _mu;

    //! Thermal conductivity
    std::unique_ptr<Conductivity> _conductivity;

    std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > _diffusivity;

    /* Below we will specialize the specialized_build_* functions to the appropriate type.
       This way, we can control how the cached transport objects get constructed
       based on the template type. This is achieved by the dummy types forcing operator
       overloading for each of the specialized types. */
    void build_conductivity( const GetPot & input, const std::string & material )
    { specialized_build_conductivity( input, material, _conductivity, conductivity_type<Conductivity>() ); }

  private:

    AntiochConstantTransportMixture();

    void specialized_build_conductivity( const GetPot & input, const std::string & material,
                                         std::unique_ptr<ConstantConductivity> & conductivity,
                                         conductivity_type<ConstantConductivity> )
    {
      conductivity.reset( new ConstantConductivity(input,material) );
    }

    void specialized_build_conductivity( const GetPot & input, const std::string & material,
                                         std::unique_ptr<ConstantPrandtlConductivity> & conductivity,
                                         conductivity_type<ConstantPrandtlConductivity> )
    {
      conductivity.reset( new ConstantPrandtlConductivity(input,material) );
    }

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename KineticsThermoCurveFit, typename Conductivity>
  inline
  libMesh::Real
  AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>::mu() const
  {
    return (*_mu)();
  }

  template<typename KineticsThermoCurveFit, typename Conductivity>
  inline
  const Conductivity&
  AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>::conductivity() const
  {
    return *_conductivity.get();
  }

  template<typename KineticsThermoCurveFit, typename Conductivity>
  inline
  const Antioch::ConstantLewisDiffusivity<libMesh::Real>&
  AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>::diffusivity() const
  {
    return *_diffusivity.get();
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_H
