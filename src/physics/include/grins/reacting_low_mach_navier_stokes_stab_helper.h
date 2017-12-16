//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef REACTING_LOW_MACH_NAVIER_STOKES_STAB_HELPER_H
#define REACTING_LOW_MACH_NAVIER_STOKES_STAB_HELPER_H

// GRINS
#include "grins/low_mach_navier_stokes_stab_helper.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesStabilizationHelper : public LowMachNavierStokesStabilizationHelper
  {
  public:

    ReactingLowMachNavierStokesStabilizationHelper( const std::string& helper_name, const GetPot& input )
      : LowMachNavierStokesStabilizationHelper(helper_name,input)
    {}

    ~ReactingLowMachNavierStokesStabilizationHelper(){};

    libMesh::Real compute_tau_species( AssemblyContext& c,
                                       unsigned int qp,
                                       libMesh::RealGradient& g,
                                       libMesh::RealTensor& G,
                                       libMesh::Real rho,
                                       libMesh::Gradient U,
                                       libMesh::Real D_s,
                                       bool is_steady ) const;

  }; // class ReactingLowMachNavierStokesStabilizationHelper

  /* ------------- Inline Functions ---------------*/
  inline
  libMesh::Real ReactingLowMachNavierStokesStabilizationHelper::compute_tau_species( AssemblyContext& c,
                                                                                     unsigned int qp,
                                                                                     libMesh::RealGradient& g,
                                                                                     libMesh::RealTensor& G,
                                                                                     libMesh::Real rho,
                                                                                     libMesh::Gradient U,
                                                                                     libMesh::Real D_s,
                                                                                     bool is_steady ) const
  {
    /*
      libMesh::Real tau = (rho*U)*(G*(rho*U)) + this->_C*D_s*D_s*G.contract(G);

      if(!is_steady)
      tau += (2.0*rho/c.get_deltat_value())*(2.0*rho/c.get_deltat_value());

      return this->_tau_factor/std::sqrt(tau);
    */
    return this->compute_tau( c, qp, D_s*D_s, g, G, rho, U, is_steady );
  }

}
#endif // REACTING_LOW_MACH_NAVIER_STOKES_STAB_HELPER_H
