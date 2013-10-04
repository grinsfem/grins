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


#ifndef GRINS_BOUSSINESQ_BUOYANCY_SPGSM_STAB_H
#define GRINS_BOUSSINESQ_BUOYANCY_SPGSM_STAB_H

// GRINS
#include "grins/boussinesq_buoyancy_base.h"
#include "grins/heat_transfer_stabilization_base.h"

namespace GRINS
{  
  //! Adds Boussinesq bouyancy adjoint stabilization source term
  /*!
    This class implements the adjiont stabilization term for the BoussinesqBuoyancy
    Physics. Intended to be used with IncompressibleNavierStokesSPGSMStabilization
    and HeatTransferSPGSMStabilization.
   */
  class BoussinesqBuoyancySPGSMStabilization : public BoussinesqBuoyancyBase,
                                               public HeatTransferStabilizationBase,
                                               public IncompressibleNavierStokesStabilizationBase
  {
  public:
    
    BoussinesqBuoyancySPGSMStabilization( const std::string& physics_name, const GetPot& input );

    ~BoussinesqBuoyancySPGSMStabilization();

    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

  protected:
    
    std::string _p_var_name;

    VariableIndex _p_var; /* Index for pressure field */

    libMeshEnums::FEFamily _P_FE_family;

    libMeshEnums::Order _P_order;

  private:

    BoussinesqBuoyancySPGSMStabilization();

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_SPGSM_STAB_H
