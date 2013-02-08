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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

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

    ~LowMachNavierStokes();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

    virtual void side_time_derivative( bool compute_jacobian,
				       libMesh::FEMContext& context,
				       CachedValues& cache );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context,
				CachedValues& cache );

    virtual void compute_element_time_derivative_cache( const libMesh::FEMContext& context,
							CachedValues& cache ) const;

    virtual void compute_element_cache( const libMesh::FEMContext& context,
					const std::vector<libMesh::Point>& points,
					CachedValues& cache ) const;

  protected:

    //! Enable pressure pinning
    bool _pin_pressure;
    
    PressurePinning _p_pinning;

    //! Helper function
    void assemble_mass_time_deriv( bool compute_jacobian, 
				   libMesh::FEMContext& context,
				   CachedValues& cache );

    //! Helper function
    void assemble_momentum_time_deriv( bool compute_jacobian, 
				       libMesh::FEMContext& context,
				       CachedValues& cache );

    //! Helper function
    void assemble_energy_time_deriv( bool compute_jacobian, 
				     libMesh::FEMContext& context,
				     CachedValues& cache );

    //! Helper function
    void assemble_continuity_mass_residual( bool compute_jacobian, 
					    libMesh::FEMContext& c );

    //! Helper function
    void assemble_momentum_mass_residual( bool compute_jacobian, 
					  libMesh::FEMContext& c );

    //! Helper function
    void assemble_energy_mass_residual( bool compute_jacobian, 
					libMesh::FEMContext& c );

    void assemble_thermo_press_elem_time_deriv( bool compute_jacobian, 
						libMesh::FEMContext& c );

    void assemble_thermo_press_side_time_deriv( bool compute_jacobian, 
						libMesh::FEMContext& c );

    void assemble_thermo_press_mass_residual( bool compute_jacobian, 
					      libMesh::FEMContext& c );

  private:

    LowMachNavierStokes();

  };

} //End namespace block

#endif // GRINS_LOW_MACH_NAVIER_STOKES_H
