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


#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_ABSTRACT_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_ABSTRACT_H

// GRINS
#include "grins_config.h"
#include "grins/grins_enums.h"
#include "grins/physics.h"
#include "grins/pressure_pinning.h"
#include "grins/assembly_context.h"
#include "grins/multi_component_vector_variable.h"

#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesAbstract : public Physics
  {
  public:

    ReactingLowMachNavierStokesAbstract(const PhysicsName& physics_name, const GetPot& input);

    virtual ~ReactingLowMachNavierStokesAbstract(){};

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    unsigned int n_species() const;

    libMesh::Real T( const libMesh::Point& p, const AssemblyContext& c ) const;

    void mass_fractions( const libMesh::Point& p, const AssemblyContext& c,
                         std::vector<libMesh::Real>& mass_fracs ) const;

    libMesh::Real rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const;

    libMesh::Real get_p0_steady( const AssemblyContext& c, unsigned int qp ) const;

    libMesh::Real get_p0_steady_side( const AssemblyContext& c, unsigned int qp ) const;

    libMesh::Real get_p0_steady( const AssemblyContext& c, const libMesh::Point& p ) const;

    libMesh::Real get_p0_transient( const AssemblyContext& c, unsigned int qp ) const;

  protected:

    libMesh::Number _p0;

    VelocityVariable& _flow_vars;
    PressureFEVariable& _press_var;
    PrimitiveTempFEVariables& _temp_vars;

    SpeciesMassFractionsVariable& _species_vars;

    /*! \todo When we mandate C++11, switch this to a std::shared_ptr. Then, in the VariableWarhouse,
      we can use dynamic_pointer_cast to get a std::shared_ptr. */
    ThermoPressureVariable* _p0_var;

    //! Number of species
    unsigned int _n_species;

    //! Gravity vector
    libMesh::Point _g;

    //! Flag to enable thermodynamic pressure calculation
    bool _enable_thermo_press_calc;

    bool _fixed_density;

    libMesh::Real _fixed_rho_value;

  private:

    ReactingLowMachNavierStokesAbstract();

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );

  }; // class ReactingLowMachNavierStokesAbstract


  inline
  unsigned int ReactingLowMachNavierStokesAbstract::n_species() const
  { return _n_species; }

  inline
  libMesh::Real ReactingLowMachNavierStokesAbstract::T( const libMesh::Point& p,
                                                        const AssemblyContext& c ) const
  { return c.point_value(_temp_vars.T(),p); }

  inline
  void ReactingLowMachNavierStokesAbstract::mass_fractions( const libMesh::Point& p,
                                                            const AssemblyContext& c,
                                                            std::vector<libMesh::Real>& mass_fracs ) const
  {
    libmesh_assert_equal_to(mass_fracs.size(), this->_n_species);

    for( unsigned int var = 0; var < this->_n_species; var++ )
      {
        mass_fracs[var] = c.point_value(_species_vars.species(var),p);
      }
  }

  inline
  libMesh::Real ReactingLowMachNavierStokesAbstract::rho( libMesh::Real T,
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
  libMesh::Real ReactingLowMachNavierStokesAbstract::get_p0_steady( const AssemblyContext& c,
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

  inline
  libMesh::Real ReactingLowMachNavierStokesAbstract::get_p0_steady_side( const AssemblyContext& c,
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

  inline
  libMesh::Real ReactingLowMachNavierStokesAbstract::get_p0_steady( const AssemblyContext& c,
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

  inline
  libMesh::Real ReactingLowMachNavierStokesAbstract::get_p0_transient( const AssemblyContext& c,
                                                                       unsigned int qp ) const
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

#endif //GRINS_REACTING_LOW_MACH_NAVIER_STOKES_ABSTRACT_H
