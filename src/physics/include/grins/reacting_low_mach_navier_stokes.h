//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_H

// GRINS
#include "grins/reacting_low_mach_navier_stokes_base.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  class ReactingLowMachNavierStokes : public ReactingLowMachNavierStokesBase<Mixture,Evaluator>
  {
  public:

    ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input);
    ~ReactingLowMachNavierStokes();
    
    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Register postprocessing variables for ReactingLowMachNavierStokes
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

    virtual void side_time_derivative( bool compute_jacobian,
				       AssemblyContext& context,
				       CachedValues& cache );

    virtual void side_constraint( bool compute_jacobian,
                                  AssemblyContext& context,
                                  CachedValues& cache );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				AssemblyContext& context,
				CachedValues& cache );

    virtual void compute_element_time_derivative_cache( const AssemblyContext& context, 
							CachedValues& cache );

    virtual void compute_side_time_derivative_cache( const AssemblyContext& context, 
						     CachedValues& cache );

    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );



  protected:

    //! Enable pressure pinning
    bool _pin_pressure;
    
    PressurePinning _p_pinning;

    //! Index from registering this quantity
    unsigned int _rho_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _species_viscosity;

    //! Index from registering this quantity
    unsigned int _mu_index;

    //! Index from registering this quantity
    unsigned int _k_index;

    //! Index from registering this quantity
    unsigned int _cp_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _mole_fractions_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _h_s_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _omega_dot_index;

  private:

    ReactingLowMachNavierStokes();

  };

} // namespace GRINS

#endif //GRINS_REACTING_LOW_MACH_NAVIER_STOKES_H
