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

#ifndef GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
#define GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_evaluator.h"
#include "grins/antioch_wilke_transport_mixture.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/wilke_evaluator.h"

namespace GRINS
{
  template<typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  class AntiochWilkeTransportEvaluator : public AntiochEvaluator<Thermo>
  {
  public:
    
    AntiochWilkeTransportEvaluator( const AntiochWilkeTransportMixture<Thermo,Viscosity,Conductivity,Diffusivity>& mixture );

    virtual ~AntiochWilkeTransportEvaluator();
    
    libMesh::Real mu( const CachedValues& cache, unsigned int qp );

    libMesh::Real k( const CachedValues& cache, unsigned int qp );

    void mu_and_k( const CachedValues& cache, unsigned int qp,
                   libMesh::Real& mu, libMesh::Real& k );

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D );

    libMesh::Real mu( const libMesh::Real T,
                      const std::vector<libMesh::Real>& Y );

    libMesh::Real k( const libMesh::Real T,
                     const std::vector<libMesh::Real>& Y );

    void D( const libMesh::Real rho, const libMesh::Real cp,
            const libMesh::Real k,
	    std::vector<libMesh::Real>& D );
    
  protected:

    boost::scoped_ptr<Antioch::WilkeEvaluator<Viscosity,Conductivity> > _wilke_evaluator;

    const Diffusivity& _diffusivity;

  private:

    AntiochWilkeTransportEvaluator();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
