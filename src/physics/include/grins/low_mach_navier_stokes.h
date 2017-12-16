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


#ifndef GRINS_LOW_MACH_NAVIER_STOKES_H
#define GRINS_LOW_MACH_NAVIER_STOKES_H

// GRINS
#include "grins/low_mach_navier_stokes_base.h"
#include "grins/pressure_pinning.h"

namespace GRINS
{

  //! Physics class for Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Incompressible Navier-Stokes equations.
  */
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokes : public LowMachNavierStokesBase<Viscosity,SpecificHeat,ThermalConductivity>
  {
  public:

    LowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input);

    ~LowMachNavierStokes(){};

    virtual void auxiliary_init( MultiphysicsSystem& system );

    //! Register postprocessing variables for LowMachNavierStokes
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    // Context initialization
    virtual void init_context( AssemblyContext & context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    virtual void compute_element_time_derivative_cache( AssemblyContext & context );

    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

  protected:

    //! Enable pressure pinning
    bool _pin_pressure;

    PressurePinning _p_pinning;

    //! Cache index for density post-processing
    unsigned int _rho_index;

    //! Helper function
    void assemble_mass_time_deriv( bool compute_jacobian,
                                   AssemblyContext& context,
                                   const CachedValues & cache );

    //! Helper function
    void assemble_momentum_time_deriv( bool compute_jacobian,
                                       AssemblyContext& context,
                                       const CachedValues & cache );

    //! Helper function
    void assemble_energy_time_deriv( bool compute_jacobian,
                                     AssemblyContext& context,
                                     const CachedValues & cache );

    //! Helper function
    void assemble_continuity_mass_residual( bool compute_jacobian,
                                            AssemblyContext & context );

    //! Helper function
    void assemble_momentum_mass_residual( bool compute_jacobian,
                                          AssemblyContext & context );

    //! Helper function
    void assemble_energy_mass_residual( bool compute_jacobian,
                                        AssemblyContext & context );

    void assemble_thermo_press_elem_time_deriv( bool compute_jacobian,
                                                AssemblyContext & context );

    void assemble_thermo_press_side_time_deriv( bool compute_jacobian,
                                                AssemblyContext & context );

    void assemble_thermo_press_mass_residual( bool compute_jacobian,
                                              AssemblyContext & context );

  private:

    LowMachNavierStokes();

  };

} //End namespace block

#endif // GRINS_LOW_MACH_NAVIER_STOKES_H
