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


#include "grins_config.h"

// This class
#include "grins/reacting_low_mach_navier_stokes_abstract.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/cantera_mixture.h"
#include "grins/grins_enums.h"
#include "grins/antioch_mixture.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  ReactingLowMachNavierStokesAbstract::ReactingLowMachNavierStokesAbstract(const std::string& physics_name,
                                                                           const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _species_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _p0_var(NULL),
      _n_species(_species_vars.n_species()),
      _fixed_density( input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/fixed_density", false ) ),
      _fixed_rho_value(0.0)
  {
    _press_var.set_is_constraint_var(true);

    this->set_parameter
      (_fixed_rho_value, input,
       "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/fixed_rho_value", 0.0 );

    _enable_thermo_press_calc = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/enable_thermo_press_calc", false );
    if( _enable_thermo_press_calc )
      {
        _p0_var = &GRINSPrivate::VariableWarehouse::get_variable_subclass<ThermoPressureVariable>(VariablesParsing::thermo_press_variable_name(input,physics_name,VariablesParsing::PHYSICS));
        _p0_var->set_is_constraint_var(true);
      }

    this->read_input_options(input);

    this->check_var_subdomain_consistency(_flow_vars);
    this->check_var_subdomain_consistency(_press_var);
    this->check_var_subdomain_consistency(_temp_vars);
    this->check_var_subdomain_consistency(_species_vars);
  }

  void ReactingLowMachNavierStokesAbstract::read_input_options( const GetPot& input )
  {
    // Read thermodynamic pressure info
    MaterialsParsing::read_property( input,
                                     "ThermodynamicPressure",
                                     PhysicsNaming::reacting_low_mach_navier_stokes(),
                                     (*this),
                                     _p0 );

    // Read gravity vector
    unsigned int g_dim = input.vector_variable_size("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g");

    _g(0) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 0 );

    if( g_dim > 1)
      _g(1) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 1 );

    if( g_dim == 3)
      _g(2) = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/g", 0.0, 2 );

  }

  void ReactingLowMachNavierStokesAbstract::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	system->time_evolving( _species_vars.species(i), 1 );
      }

    system->time_evolving(_flow_vars.u(),1);

    if (dim > 1)
      system->time_evolving(_flow_vars.v(),1);

    if (dim == 3)
      system->time_evolving(_flow_vars.w(),1);

    system->time_evolving(_temp_vars.T(),1);
    system->time_evolving(_press_var.p(),1);

    if( _enable_thermo_press_calc )
      system->time_evolving(_p0_var->p0(),1);
  }

  void ReactingLowMachNavierStokesAbstract::init_context( AssemblyContext& context )
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
  }

} // end namespace GRINS
