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
#include "grins/antioch_constant_transport_evaluator.h"

// GRINS
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename KT, typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<KT,Thermo,Conductivity>::mu( const libMesh::Real /*T*/,
                                                                               const libMesh::Real /*P*/,
                                                                               const std::vector<libMesh::Real>& /*Y*/ )
  {
    return _mu;
  }

  template<typename KT, typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<KT,Thermo,Conductivity>::k( const libMesh::Real T,
                                                                              const libMesh::Real /*P*/,
                                                                              const std::vector<libMesh::Real>& Y )
  {
    // Second T is dummy
    const libMesh::Real cp = this->cp( T, T, Y );

    return _conductivity( _mu, cp );
  }

  template<typename KT, typename Thermo, typename Conductivity>
  void AntiochConstantTransportEvaluator<KT,Thermo,Conductivity>::mu_and_k_and_D( const libMesh::Real /*T*/,
                                                                                  const libMesh::Real rho,
                                                                                  const libMesh::Real cp,
                                                                                  const std::vector<libMesh::Real>& /*Y*/,
                                                                                  libMesh::Real& mu, libMesh::Real& k,
                                                                                  std::vector<libMesh::Real>& D )
  {
    mu = _mu;

    k = _conductivity( _mu, cp );

    std::fill( D.begin(), D.end(), _diffusivity.D(rho,cp,k) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
