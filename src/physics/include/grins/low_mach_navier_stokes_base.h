//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_LOW_MACH_NAVIER_STOKES_BASE_H
#define GRINS_LOW_MACH_NAVIER_STOKES_BASE_H

// GRINS
#include "grins/physics.h"
#include "grins/assembly_context.h"
#include "grins/grins_enums.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/point.h"

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
    virtual void init_context( AssemblyContext& context );

    libMesh::Real T( const libMesh::Point& p, const AssemblyContext& c ) const;

    libMesh::Real rho( libMesh::Real T, libMesh::Real p0 ) const;
    
    libMesh::Real d_rho_dT( libMesh::Real T, libMesh::Real p0 ) const;

    libMesh::Real get_p0_steady( const AssemblyContext& c, unsigned int qp ) const;

    libMesh::Real get_p0_steady_side( const AssemblyContext& c, unsigned int qp ) const;

    libMesh::Real get_p0_steady( const AssemblyContext& c, const libMesh::Point& p ) const;

    libMesh::Real get_p0_transient( AssemblyContext& c, unsigned int qp ) const;

  protected:

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
    GRINSEnums::FEFamily _V_FE_family, _P_FE_family, _T_FE_family;

    //! Element orders, read from input
    GRINSEnums::Order _V_order, _P_order, _T_order;

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

  template<class V, class SH, class TC>
  inline
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::T( const libMesh::Point& p, const AssemblyContext& c ) const
  {
    return c.point_value(_T_var,p);
  }

  template<class V, class SH, class TC>
  inline
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::rho( libMesh::Real T, libMesh::Real p0 ) const
  {
    return p0/(this->_R*T);
  }
  
  template<class V, class SH, class TC>
  inline
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::d_rho_dT( libMesh::Real T, libMesh::Real p0 ) const
  {
    return -p0/(this->_R*(T*T));
  }
  
  template<class V, class SH, class TC>
  inline 
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::get_p0_steady( const AssemblyContext& c,
								 unsigned int qp ) const
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

  template<class V, class SH, class TC>
  inline 
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::get_p0_steady_side( const AssemblyContext& c,
								      unsigned int qp ) const
  {
    libMesh::Real p0;
    if( this->_enable_thermo_press_calc )
      {
	p0 = c.side_value( _p0_var, qp );
      }
    else
      {
	p0 = _p0;
      }
    return p0;
  }

  template<class V, class SH, class TC>
  inline 
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::get_p0_steady( const AssemblyContext& c, 
								 const libMesh::Point& p ) const
  {
    libMesh::Real p0;
    if( this->_enable_thermo_press_calc )
      {
	p0 = c.point_value( _p0_var, p );
      }
    else
      {
	p0 = _p0;
      }
    return p0;
  }

  template<class V, class SH, class TC>
  inline 
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::get_p0_transient( AssemblyContext& c, unsigned int qp ) const
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

} // namespace GRINS

#endif // GRINS_LOW_MACH_NAVIER_STOKES_BASE_H
