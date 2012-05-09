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

#ifndef LOW_MACH_NAVIER_STOKES_H
#define LOW_MACH_NAVIER_STOKES_H

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
#include "config.h"
#include "physics.h"
#include "constant_viscosity.h"
#include "constant_specific_heat.h"
#include "constant_conductivity.h"

namespace GRINS
{

  //! Physics class for Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Incompressible Navier-Stokes equations.
   */
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokes : public Physics
  {
  public:

    LowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input);

    ~LowMachNavierStokes();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization of Navier-Stokes variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::DiffContext &context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    // Constraint part(s)
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );
    //! Handles Dirichlet boundary conditions
    /*! Note that for any generic function specifications, 
      any components not specified will be assigned a zero Dirichlet value. */
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    // Mass matrix part(s)
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

    virtual int string_to_int( const std::string& bc_type_in );

    virtual void init_bc_data( const GRINS::BoundaryID bc_id, 
			       const std::string& bc_id_string, 
			       const int bc_type, 
			       const GetPot& input );

    virtual void init_dirichlet_bcs( libMesh::DofMap& dof_map );

  protected:

    //! Thermodynamic pressure divided by gas constant
    libMesh::Number _p0_over_R;

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

    //! Used for storing values corresponding to GRINS::PRESCRIBED_VELOCITY values
    std::map< unsigned int, std::vector<double> > _vel_boundary_values;

    //! Enable pressure pinning
    bool _pin_pressure;
    
    //! Value of pressure we wish to pin
    libMesh::Number _pin_value;

    //! Location we want to pin the pressure
    libMesh::Point _pin_location;

    enum LMNS_BC_TYPES{NO_SLIP=0, PRESCRIBED_VELOCITY, INFLOW, ISOTHERMAL_WALL, ADIABATIC_WALL, PRESCRIBED_HEAT_FLUX, GENERAL_HEAT_FLUX};

    //! Helper function
    void assemble_mass_time_deriv( bool request_jacobian, 
				   libMesh::FEMContext& context, 
				   libMesh::FEMSystem* system );

    //! Helper function
    void assemble_momentum_time_deriv( bool request_jacobian, 
				       libMesh::FEMContext& context, 
				       libMesh::FEMSystem* system );

    //! Helper function
    void assemble_energy_time_deriv( bool request_jacobian, 
				     libMesh::FEMContext& context, 
				     libMesh::FEMSystem* system );

    //! Helper function
    void assemble_continuity_mass_residual( bool request_jacobian, 
					    libMesh::FEMContext& c, 
					    libMesh::FEMSystem* system );

    //! Helper function
    void assemble_momentum_mass_residual( bool request_jacobian, 
					  libMesh::FEMContext& c, 
					  libMesh::FEMSystem* system );

    //! Helper function
    void assemble_energy_mass_residual( bool request_jacobian, 
					libMesh::FEMContext& c, 
					libMesh::FEMSystem* system );

  private:
    LowMachNavierStokes();

  };

} //End namespace block

#endif // LOW_MACH_NAVIER_STOKES_H
