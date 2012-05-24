//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef LOW_MACH_NAVIER_STOKES_STAB_BASE_H
#define LOW_MACH_NAVIER_STOKES_STAB_BASE_H

//GRINS
#include "low_mach_navier_stokes_base.h"
#include "stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesStabilizationBase : public LowMachNavierStokesBase<Viscosity,SpecificHeat,ThermalConductivity>
  {

  public:

    LowMachNavierStokesStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~LowMachNavierStokesStabilizationBase();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( const GetPot& input );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::DiffContext &context );

    inline
    libMesh::Real compute_tau_continuity( libMesh::Real tau_M,
					  libMesh::RealGradient& g,
					  libMesh::RealTensor&,
					  libMesh::Gradient,
					  libMesh::Real  ) const
    {
      return this->_tau_factor/(tau_M*(g*g));
    }

    inline
    libMesh::Real compute_tau_momentum( libMesh::FEMContext& c,
					unsigned int qp,
					libMesh::RealGradient& g,
					libMesh::RealTensor& G,
					libMesh::Real rho,
					libMesh::Gradient U,
					libMesh::Real T,
					bool is_steady ) const
    {
      libMesh::Real mu = this->_mu(T);

      return this->compute_tau( c, qp, mu*mu, g, G, rho, U, is_steady );
    }

    inline
    libMesh::Real compute_tau_energy( libMesh::FEMContext& c,
				      unsigned int qp,
				      libMesh::RealGradient& g,
				      libMesh::RealTensor& G,
				      libMesh::Real rho,
				      libMesh::Gradient U,
				      libMesh::Real T,
				      bool is_steady ) const
    {
      libMesh::Real k = this->_k(T);
      libMesh::Real cp = this->_cp(T);

      return this->compute_tau( c, qp, k*k, g, G, rho*cp, U, is_steady );
    }

    inline
    libMesh::Real compute_tau( libMesh::FEMContext& c,
			       unsigned int qp,
			       libMesh::Real mat_prop_sq,
			       libMesh::RealGradient& g,
			       libMesh::RealTensor& G,
			       libMesh::Real rho,
			       libMesh::Gradient U,
			       bool is_steady ) const
    {
      libMesh::Real tau = (rho*U)*(G*(rho*U)) + this->_C*mat_prop_sq*G.contract(G);

      if(!is_steady)
	tau += (2.0*rho/c.get_deltat_value())*(2.0*rho/c.get_deltat_value());

      return this->_tau_factor/std::sqrt(tau);
    }

    libMesh::Real compute_res_continuity_steady( libMesh::FEMContext& context,
						 unsigned int qp ) const;
    
    libMesh::Real compute_res_continuity_transient( libMesh::FEMContext& context,
						    unsigned int qp ) const;
    
    libMesh::RealGradient compute_res_momentum_steady( libMesh::FEMContext& context,
						       unsigned int qp ) const;

    libMesh::RealGradient compute_res_momentum_transient( libMesh::FEMContext& context,
							  unsigned int qp ) const;

    libMesh::Real compute_res_energy_steady( libMesh::FEMContext& context,
					     unsigned int qp ) const;
    
    libMesh::Real compute_res_energy_transient( libMesh::FEMContext& context,
						unsigned int qp ) const;


  protected:

    libMesh::Real _C, _tau_factor;

    StabilizationHelper _stab_helper;
    
  private:
    LowMachNavierStokesStabilizationBase();

  }; // End LowMachNavierStokesStabilizationBase class declarations

} // End namespace GRINS

#endif //LOW_MACH_NAVIER_STOKES_STAB_BASE_H
