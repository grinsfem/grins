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
#include "grins/antioch_constant_transport_mixture.h"

// GRINS
#include "grins/materials_parsing.h"

namespace GRINS
{
  template<typename KineticsThermoCurveFit,typename Conductivity>
  AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>::
  AntiochConstantTransportMixture( const GetPot& input,const std::string& material )
    : AntiochMixture<KineticsThermoCurveFit>(input,material)
  {
    {
      std::string warning = "==============================================\n";
      warning += "WARNING: This AntiochConstantTransportMixture constructor is DEPREACTED!\n";
      warning += "         Prefer alternate constructor where parsing\n";
      warning += "         is done outside this class.\n";
      warning += "==============================================\n";

      libmesh_warning(warning);
    }

    libMesh::Real Le = MaterialsParsing::parse_lewis_number(input,material);
    _diffusivity.reset( new Antioch::ConstantLewisDiffusivity<libMesh::Real>(Le) );
    _mu.reset( new ConstantViscosity(input,material) );
    this->build_conductivity(input,material);
  }

  template<typename KineticsThermoCurveFit,typename Conductivity>
  AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>::
  AntiochConstantTransportMixture
  ( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > & nasa_mixture,
    std::unique_ptr<ConstantViscosity> & visc,
    std::unique_ptr<Conductivity> & cond,
    std::unique_ptr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > & diff,
    libMesh::Real min_T,
    bool clip_negative_rho )
    : AntiochMixture<KineticsThermoCurveFit>(chem_mixture,reaction_set,nasa_mixture,min_T,clip_negative_rho)
  {
    /*! \todo Use std::move() when we have C++11 */
    _mu.reset( visc.release() );
    _conductivity.reset( cond.release() );
    _diffusivity.reset( diff.release() );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
