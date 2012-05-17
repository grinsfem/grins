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
#ifndef LOW_MACH_NAVIER_STOKES_VMS_STAB_H
#define LOW_MACH_NAVIER_STOKES_VMS_STAB_H

//libMesh
#include "time_solver.h"

//GRINS
#include "low_mach_navier_stokes_base.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesVMSStabilization : public LowMachNavierStokesBase<Viscosity,SpecificHeat,ThermalConductivity>
  {

  public:

    LowMachNavierStokesVMSStabilization( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~LowMachNavierStokesVMSStabilization();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( const GetPot& input );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::DiffContext &context );

    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

  protected:

    libMesh::Real _C;

    void assemble_continuity_time_deriv( bool request_jacobian,
					 libMesh::FEMContext& context,
					 libMesh::FEMSystem* system );

    void assemble_momentum_time_deriv( bool request_jacobian,
				       libMesh::FEMContext& context,
				       libMesh::FEMSystem* system );

    void assemble_energy_time_deriv( bool request_jacobian,
				     libMesh::FEMContext& context,
				     libMesh::FEMSystem* system );

    void assemble_continuity_mass_residual( bool request_jacobian,
					    libMesh::FEMContext& context,
					    libMesh::FEMSystem* system );

    void assemble_momentum_mass_residual( bool request_jacobian,
					  libMesh::FEMContext& context,
					  libMesh::FEMSystem* system );

    void assemble_energy_mass_residual( bool request_jacobian,
					libMesh::FEMContext& context,
					libMesh::FEMSystem* system );

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

    inline
    libMesh::Real compute_tau_continuity( libMesh::Real tau_M,
					  libMesh::RealGradient& g ) const
    {
      return 1.0/(tau_M*(g*g));
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
      libMesh::Real k_over_cp = (this->_k(T))/(this->_cp(T));

      return this->compute_tau( c, qp, k_over_cp*k_over_cp, g, G, rho, U, is_steady );
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

      return 1.0/std::sqrt(tau);
    }

    inline
    libMesh::RealGradient compute_g( libMesh::FEMContext& c,
				     unsigned int qp ) const
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];
      return libMesh::RealGradient( fe->get_dxidx()[qp] + fe->get_detadx()[qp] + fe->get_dzetadx()[qp],
				    fe->get_dxidy()[qp] + fe->get_detady()[qp] + fe->get_dzetady()[qp],
				    fe->get_dxidz()[qp] + fe->get_detadz()[qp] + fe->get_dzetadz()[qp] );
    }
    
    inline
    libMesh::RealTensor compute_G( libMesh::FEMContext& c,
				   unsigned int qp ) const
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];
      libMesh::Real dxidx = fe->get_dxidx()[qp];
      libMesh::Real dxidy = fe->get_dxidy()[qp];
      libMesh::Real dxidz = fe->get_dxidz()[qp];

      libMesh::Real detadx = fe->get_detadx()[qp];
      libMesh::Real detady = fe->get_detady()[qp];
      libMesh::Real detadz = fe->get_detadz()[qp];

      libMesh::Real dzetadx = fe->get_dzetadx()[qp];
      libMesh::Real dzetady = fe->get_dzetady()[qp];
      libMesh::Real dzetadz = fe->get_dzetadz()[qp];
      
      return libMesh::RealTensor( dxidx*dxidx + detadx*detadx + dzetadx*dzetadx,
				  dxidx*dxidy + detadx*detady + dzetadx*dzetady,
				  dxidx*dxidz + detadx*detadz + dzetadx*dzetadz,
				  dxidy*dxidx + detady*detadx + dzetady*dzetadx,
				  dxidy*dxidy + detady*detady + dzetady*dzetady,
				  dxidy*dxidz + detady*detadz + dzetady*dzetadz,
				  dxidz*dxidx + detadz*detadx + dzetadz*dzetadx,
				  dxidz*dxidy + detadz*detady + dzetadz*dzetady,
				  dxidz*dxidz + detadz*detadz + dzetadz*dzetadz );
    }

    
  private:
    LowMachNavierStokesVMSStabilization();

  }; // End LowMachNavierStokesVMSStabilization class declarations

} // End namespace GRINS

#endif //LOW_MACH_NAVIER_STOKES_VMS_STAB_H
