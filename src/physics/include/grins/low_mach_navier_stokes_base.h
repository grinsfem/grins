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


#ifndef GRINS_LOW_MACH_NAVIER_STOKES_BASE_H
#define GRINS_LOW_MACH_NAVIER_STOKES_BASE_H

// GRINS
#include "grins/physics.h"
#include "grins/assembly_context.h"
#include "grins/grins_enums.h"
#include "grins/multi_component_vector_variable.h"

#include "grins/single_variable.h"

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

    LowMachNavierStokesBase(const PhysicsName& physics_name, const std::string& core_physics_name, const GetPot& input);

    ~LowMachNavierStokesBase(){};

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

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    //! Thermodynamic pressure divided by gas constant
    libMesh::Number _p0_over_R;

    libMesh::Number _p0, _R, _T0;

    VelocityVariable& _flow_vars;
    PressureFEVariable& _press_var;
    PrimitiveTempFEVariables& _temp_vars;

    /*! \todo When we mandate C++11, switch this to a std::shared_ptr. Then, in the VariableWarhouse,
      we can use dynamic_pointer_cast to get a std::shared_ptr. */
    ThermoPressureVariable*  _p0_var;

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

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );
  };

  template<class V, class SH, class TC>
  inline
  libMesh::Real LowMachNavierStokesBase<V,SH,TC>::T( const libMesh::Point& p, const AssemblyContext& c ) const
  {
    return c.point_value(_temp_vars.T(),p);
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
        p0 = c.interior_value( _p0_var->p0(), qp );
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
        p0 = c.side_value( _p0_var->p0(), qp );
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
        p0 = c.point_value( _p0_var->p0(), p );
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
        p0 = c.fixed_interior_value( _p0_var->p0(), qp );
      }
    else
      {
        p0 = _p0;
      }
    return p0;
  }

} // namespace GRINS

#endif // GRINS_LOW_MACH_NAVIER_STOKES_BASE_H
