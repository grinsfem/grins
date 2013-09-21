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


#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BASE_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BASE_H

// GRINS
#include "grins_config.h"
#include "grins/physics.h"
#include "grins/pressure_pinning.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesBase : public Physics
  {
  public:

    ReactingLowMachNavierStokesBase(const PhysicsName& physics_name, const GetPot& input);
    ~ReactingLowMachNavierStokesBase();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext &context );

    unsigned int n_species() const;

    libMesh::Real T( const libMesh::Point& p, const libMesh::FEMContext& c ) const;

    void mass_fractions( const libMesh::Point& p, const libMesh::FEMContext& c,
                         std::vector<libMesh::Real>& mass_fracs ) const;

    libMesh::Real rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const;

    libMesh::Real get_p0_steady( const libMesh::FEMContext& c, unsigned int qp ) const;

    libMesh::Real get_p0_steady_side( const libMesh::FEMContext& c, unsigned int qp ) const;
 
    libMesh::Real get_p0_steady( const libMesh::FEMContext& c, const libMesh::Point& p ) const;

    libMesh::Real get_p0_transient( const libMesh::FEMContext& c, unsigned int qp ) const;

    //! Method to interface to thermochemistry quantity
    /*! Intended to be called externally from Physics class, e.g. in a QoI.
        This way, we can leverage the Physics class for the thermochemistry.*/
    virtual libMesh::Real cp_mix( const libMesh::Real T,
                                  const std::vector<libMesh::Real>& Y ) =0;

    //! Method to interface to thermochemistry quantity
    /*! Intended to be called externally from Physics class, e.g. in a QoI.
        This way, we can leverage the Physics class for the thermochemistry.*/
    virtual libMesh::Real mu( const libMesh::Real T,
                              const std::vector<libMesh::Real>& Y ) =0;

    //! Method to interface to thermochemistry quantity
    /*! Intended to be called externally from Physics class, e.g. in a QoI.
        This way, we can leverage the Physics class for the thermochemistry.*/
    virtual libMesh::Real k( const libMesh::Real T,
                             const std::vector<libMesh::Real>& Y ) =0;

    //! Method to interface to thermochemistry quantity
    /*! Intended to be called externally from Physics class, e.g. in a QoI.
        This way, we can leverage the Physics class for the thermochemistry.*/
    virtual void D( const libMesh::Real rho, const libMesh::Real cp,
                    const libMesh::Real k,
                    std::vector<libMesh::Real>& D ) =0;

  protected:

    libMesh::Number _p0;

    //! Physical dimension of problem
    unsigned int _dim;

    //! Number of species
    unsigned int _n_species;

    //! Indices for each (owned) variable;
    std::vector<VariableIndex> _species_vars; /* Indicies for species densities */
    VariableIndex _u_var; /* Index for x-velocity field */
    VariableIndex _v_var; /* Index for y-velocity field */
    VariableIndex _w_var; /* Index for z-velocity field */
    VariableIndex _p_var; /* Index for pressure field */
    VariableIndex _T_var; /* Index for pressure field */
    VariableIndex _p0_var; /* Index for thermodynamic pressure */

    //! Names of each (owned) variable in the system
    std::vector<std::string> _species_var_names;
    std::string _u_var_name, _v_var_name, _w_var_name, _p_var_name, _T_var_name, _p0_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _species_FE_family, _V_FE_family, _P_FE_family, _T_FE_family;

    //! Element orders, read from input
    libMeshEnums::Order _species_order, _V_order, _P_order, _T_order;

    //! Gravity vector
    libMesh::Point _g; 

    //! Flag to enable thermodynamic pressure calculation
    bool _enable_thermo_press_calc;

    bool _fixed_density;

    libMesh::Real _fixed_rho_value;

  private:

    ReactingLowMachNavierStokesBase();

  }; // class ReactingLowMachNavierStokesBase

  inline
  unsigned int ReactingLowMachNavierStokesBase::n_species() const
  { return _n_species; }

  
  inline
  libMesh::Real ReactingLowMachNavierStokesBase::T( const libMesh::Point& p, 
                                                    const libMesh::FEMContext& c ) const
  { return c.point_value(_T_var,p); }

  inline
  void ReactingLowMachNavierStokesBase::mass_fractions( const libMesh::Point& p, 
                                                        const libMesh::FEMContext& c,
                                                        std::vector<libMesh::Real>& mass_fracs ) const
  {
    libmesh_assert_equal_to(mass_fracs.size(), this->_n_species);

    for( unsigned int var = 0; var < this->_n_species; var++ )
      {
        mass_fracs[var] = c.point_value(_species_vars[var],p);
      }

    return;
  }

  inline
  libMesh::Real ReactingLowMachNavierStokesBase::rho( libMesh::Real T, 
                                                      libMesh::Real p0,
                                                      libMesh::Real R_mix) const
  {
    libMesh::Real value = 0;
    if( this->_fixed_density )
      value = this->_fixed_rho_value;
    else
      value = p0/(R_mix*T);

    return value;
  }

  inline
  libMesh::Real ReactingLowMachNavierStokesBase::get_p0_steady( const libMesh::FEMContext& c, 
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

  inline
  libMesh::Real ReactingLowMachNavierStokesBase::get_p0_steady_side( const libMesh::FEMContext& c, 
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

  inline
  libMesh::Real ReactingLowMachNavierStokesBase::get_p0_steady( const libMesh::FEMContext& c, 
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
  
  inline
  libMesh::Real ReactingLowMachNavierStokesBase::get_p0_transient( const libMesh::FEMContext& c,
                                                                   unsigned int qp ) const
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

#endif //GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BASE_H
