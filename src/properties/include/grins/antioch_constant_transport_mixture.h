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

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename Conductivity>
  class AntiochConstantTransportMixture : public AntiochMixture
  {
  public:

    AntiochConstantTransportMixture( const GetPot& input );

    virtual ~AntiochConstantTransportMixture();

    libMesh::Real mu() const;

    const Conductivity& conductivity() const;

    const Antioch::ConstantLewisDiffusivity<libMesh::Real>& diffusivity() const;

    typedef AntiochChemistry ChemistryParent;

  protected:

    libMesh::Real _mu;

    boost::scoped_ptr<Conductivity> _conductivity;

    boost::scoped_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > _diffusivity;

    /* Below we will specialize the specialized_build_* functions to the appropriate type.
       This way, we can control how the cached transport objects get constructed
       based on the template type. This is achieved by the dummy types forcing operator
       overloading for each of the specialized types. */
    void build_conductivity( const GetPot& input )
    { specialized_build_conductivity( input, _conductivity, conductivity_type<Conductivity>() ); }

  private:

    AntiochConstantTransportMixture();

    void specialized_build_conductivity( const GetPot& input,
                                         boost::scoped_ptr<ConstantConductivity>& conductivity,
                                         conductivity_type<ConstantConductivity> )
    {
      conductivity.reset( new ConstantConductivity(input) );
      return;
    }

    void specialized_build_conductivity( const GetPot& input,
                                         boost::scoped_ptr<ConstantPrandtlConductivity>& conductivity,
                                         conductivity_type<ConstantPrandtlConductivity> )
    {
      conductivity.reset( new ConstantPrandtlConductivity(input) );
      return;
    }

  };
  
  /* ------------------------- Inline Functions -------------------------*/
  template<typename Conductivity>
  inline
  libMesh::Real AntiochConstantTransportMixture<Conductivity>::mu() const
  {
    return _mu;
  }

  template<typename Conductivity>
  inline
  const Conductivity& AntiochConstantTransportMixture<Conductivity>::conductivity() const
  {
    return *_conductivity.get();
  }

  template<typename Conductivity>
  inline
  const Antioch::ConstantLewisDiffusivity<libMesh::Real>& AntiochConstantTransportMixture<Conductivity>::diffusivity() const
  {
    return *_diffusivity.get();
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CONSTANT_TRANSPORT_MIXTURE_H
