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


#include "grins_config.h"

// This class
#include "grins/reacting_low_mach_navier_stokes_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/cantera_mixture.h"
#include "grins/grins_enums.h"
#include "grins/antioch_mixture.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokesBase<Mixture,Evaluator>::ReactingLowMachNavierStokesBase(const std::string& physics_name,
									    const GetPot& input)
    : Physics(physics_name, input),
      _gas_mixture(input,MaterialsParsing::material_name(input,PhysicsNaming::reacting_low_mach_navier_stokes())),
      _flow_vars(input, PhysicsNaming::reacting_low_mach_navier_stokes()),
      _press_var(input,PhysicsNaming::reacting_low_mach_navier_stokes()),
      _temp_vars(input, PhysicsNaming::reacting_low_mach_navier_stokes()),
      _p0_var(input, PhysicsNaming::reacting_low_mach_navier_stokes()),
      _species_vars(input, PhysicsNaming::reacting_low_mach_navier_stokes()),
      _n_species(_species_vars.n_species()),
      _fixed_density( input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/fixed_density", false ) ),
      _fixed_rho_value(0.0)
  {
    this->set_parameter
      (_fixed_rho_value, input,
       "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/fixed_rho_value", 0.0 );

    this->read_input_options(input);

    return;
  }

  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokesBase<Mixture,Evaluator>::~ReactingLowMachNavierStokesBase()
  {
    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::read_input_options( const GetPot& input )
  {
    // Read thermodynamic pressure info
    MaterialsParsing::read_property( input,
                                     "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/p0",
                                     "ThermodynamicPressure",
                                     PhysicsNaming::reacting_low_mach_navier_stokes(),
                                     (*this),
                                     _p0 );

    _enable_thermo_press_calc = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/enable_thermo_press_calc", false );

    // Read gravity vector
    unsigned int g_dim = input.vector_variable_size("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g");

    _g(0) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 1 );

    if( g_dim == 3)
      _g(2) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 2 );

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    this->_species_vars.init(system);
    this->_flow_vars.init(system);
    this->_press_var.init(system);
    this->_temp_vars.init(system);

    /* If we need to compute the thermodynamic pressure, we force this to be a first
       order scalar variable. */
    if( _enable_thermo_press_calc )
      _p0_var.init(system);

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	system->time_evolving( _species_vars.species(i) );
      }

    system->time_evolving(_flow_vars.u());
    system->time_evolving(_flow_vars.v());

    if (dim == 3)
      system->time_evolving(_flow_vars.w());

    system->time_evolving(_temp_vars.T());
    system->time_evolving(_press_var.p());

    if( _enable_thermo_press_calc )
      system->time_evolving(_p0_var.p0());

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_species_vars.species(0))->get_JxW();
    context.get_element_fe(_species_vars.species(0))->get_phi();
    context.get_element_fe(_species_vars.species(0))->get_dphi();
    context.get_element_fe(_species_vars.species(0))->get_xyz();

    context.get_element_fe(_flow_vars.u())->get_JxW();
    context.get_element_fe(_flow_vars.u())->get_phi();
    context.get_element_fe(_flow_vars.u())->get_dphi();
    context.get_element_fe(_flow_vars.u())->get_xyz();

    context.get_element_fe(_temp_vars.T())->get_JxW();
    context.get_element_fe(_temp_vars.T())->get_phi();
    context.get_element_fe(_temp_vars.T())->get_dphi();
    context.get_element_fe(_temp_vars.T())->get_xyz();

    context.get_element_fe(_press_var.p())->get_phi();
    context.get_element_fe(_press_var.p())->get_xyz();

    return;
  }

} // end namespace GRINS
