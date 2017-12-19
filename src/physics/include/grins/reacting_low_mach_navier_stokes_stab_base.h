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

#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_STAB_BASE_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_STAB_BASE_H

//GRINS
#include "grins/reacting_low_mach_navier_stokes_base.h"
#include "grins/reacting_low_mach_navier_stokes_stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<typename Mixture, typename Evaluator>
  class ReactingLowMachNavierStokesStabilizationBase : public ReactingLowMachNavierStokesBase<Mixture>
  {

  public:

    ReactingLowMachNavierStokesStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                                  std::unique_ptr<Mixture> & gas_mix);

    virtual ~ReactingLowMachNavierStokesStabilizationBase(){};

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    void compute_res_steady( AssemblyContext& context,
                             unsigned int qp,
                             libMesh::Real& RP_s,
                             libMesh::RealGradient& RM_s,
                             libMesh::Real& RE_s,
                             std::vector<libMesh::Real>& Rs_s );

    void compute_res_transient( AssemblyContext& context,
                                unsigned int qp,
                                libMesh::Real& RP_t,
                                libMesh::RealGradient& RM_t,
                                libMesh::Real& RE_t,
                                std::vector<libMesh::Real>& Rs_t );

  protected:

    ReactingLowMachNavierStokesStabilizationHelper _stab_helper;

  private:

    ReactingLowMachNavierStokesStabilizationBase();

  };

} // end namespace GRINS

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_STAB_BASE_H
