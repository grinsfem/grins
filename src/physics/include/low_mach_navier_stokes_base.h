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

#ifndef LOW_MACH_NAVIER_STOKES_BASE_H
#define LOW_MACH_NAVIER_STOKES_BASE_H

// libMesh
#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"
#include "string_to_enum.h"
#include "fem_system.h"
#include "fem_context.h"

// GRINS
#include "grins_config.h"
#include "physics.h"
#include "pressure_pinning.h"
#include "constant_viscosity.h"
#include "constant_specific_heat.h"
#include "constant_conductivity.h"
#include "grins/low_mach_navier_stokes_bc_handling.h"

namespace GRINS
{

  //! Physics class for Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Incompressible Navier-Stokes equations.
   */
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesBase : public Physics
  {
  public:

    LowMachNavierStokesBase(const PhysicsName& physics_name, const GetPot& input);

    ~LowMachNavierStokesBase();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

  protected:

    inline
    libMesh::Real compute_rho( libMesh::Real T, libMesh::Real p0 ) const
    {
      return p0/(this->_R*T);
    }

    inline 
    libMesh::Real get_p0_steady( libMesh::FEMContext& c, unsigned int qp ) const
    {
      libMesh::Real p0;
      if( this->_enable_thermo_press_calc )
	{
	  p0 = c.interior_value( _p0_var, qp );
	}
      else
	{
	  p0 = _p0;
	}
      return p0;
    }

    inline 
    libMesh::Real get_p0_transient( libMesh::FEMContext& c, unsigned int qp ) const
    {
      libMesh::Real p0;
      if( this->_enable_thermo_press_calc )
	{
	  p0 = c.fixed_interior_value( _p0_var, qp );
	}
      else
	{
	  p0 = _p0;
	}
      return p0;
    }

    //! Thermodynamic pressure divided by gas constant
    libMesh::Number _p0_over_R;

    libMesh::Number _p0, _R, _T0;

    //! Physical dimension of problem
    unsigned int _dim;

    //! Indices for each (owned) variable;
    VariableIndex _u_var; /* Index for x-velocity field */
    VariableIndex _v_var; /* Index for y-velocity field */
    VariableIndex _w_var; /* Index for z-velocity field */
    VariableIndex _p_var; /* Index for pressure field */
    VariableIndex _T_var; /* Index for pressure field */
    VariableIndex _p0_var; /* Index for thermodynamic pressure */

    //! Names of each (owned) variable in the system
    std::string _u_var_name, _v_var_name, _w_var_name, _p_var_name, _T_var_name, _p0_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _V_FE_family, _P_FE_family, _T_FE_family;

    //! Element orders, read from input
    libMeshEnums::Order _V_order, _P_order, _T_order;

    //! Viscosity object
    Viscosity _mu;

    //! Specific heat object
    SpecificHeat _cp;

    //! Thermal conductivity object
    ThermalConductivity _k;

    //! Gravity vector
    libMesh::Point _g; 

    //! Flag to enable thermodynamic pressure calculation
    bool _enable_thermo_press_calc;

  private:
    LowMachNavierStokesBase();

  };

} //End namespace block

#endif // LOW_MACH_NAVIER_STOKES_BASE_H
