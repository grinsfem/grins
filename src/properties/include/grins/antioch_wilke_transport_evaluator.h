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

namespace GRINS
{
  template<typename Viscosity, typename Conductivity>
  class AntiochWilkeTransportEvaluator
  {
  public:
    
    AntiochWilkeTransportEvaluator( const AntiochMixture& mixture );
    ~AntiochWilkeTransportEvaluator();

    libMesh::Real mu( const libMesh::Real T, const std::vector<libMesh::Real>& mass_fractions );

    libMesh::Real k( const libMesh::Real T, const std::vector<libMesh::Real>& mass_fractions );

    void mu_and_k( const libMesh::Real T, const std::vector<libMesh::Real>& mass_fractions,
                   libMesh::Real& mu, libMesh::Real& k );

  protected:

    Antioch::WilkeEvaluator<Viscosity,Conductivity> _wilke_evaluator;

  private:

    AntiochWilkeTransportEvaluator();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
