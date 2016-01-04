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
      _gas_mixture(input,MaterialsParsing::material_name(input,reacting_low_mach_navier_stokes)),
      _flow_vars(input, reacting_low_mach_navier_stokes),
      _temp_vars(input, reacting_low_mach_navier_stokes),
      _p0_var(input, reacting_low_mach_navier_stokes),
      _fixed_density( input("Physics/"+reacting_low_mach_navier_stokes+"/fixed_density", false ) ),
      _fixed_rho_value(0.0)
  {
    this->set_parameter
      (_fixed_rho_value, input,
       "Physics/"+reacting_low_mach_navier_stokes+"/fixed_rho_value", 0.0 );

    // Parse species and setup variable names
    MaterialsParsing::parse_species_varnames(input,MaterialsParsing::material_name(input,reacting_low_mach_navier_stokes),_species_var_names);
    this->_n_species = _species_var_names.size();

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
    // Read FE family info
    this->_species_FE_family = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_FE_family", "LAGRANGE") );

    // Read FE family info
    this->_species_order = libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_order", "SECOND") );

    // Read thermodynamic pressure info
    MaterialsParsing::read_property( input,
                                     "Physics/"+reacting_low_mach_navier_stokes+"/p0",
                                     "ThermodynamicPressure",
                                     reacting_low_mach_navier_stokes,
                                     (*this),
                                     _p0 );

    _enable_thermo_press_calc = input("Physics/"+reacting_low_mach_navier_stokes+"/enable_thermo_press_calc", false );

    // Read gravity vector
    unsigned int g_dim = input.vector_variable_size("Physics/"+reacting_low_mach_navier_stokes+"/g");

    _g(0) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 1 );

    if( g_dim == 3)
      _g(2) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 2 );

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _species_vars.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	_species_vars.push_back( system->add_variable( _species_var_names[i],
						       this->_species_order, _species_FE_family) );
      }

    this->_flow_vars.init(system);
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
	system->time_evolving( _species_vars[i] );
      }

    system->time_evolving(_flow_vars.u_var());
    system->time_evolving(_flow_vars.v_var());

    if (dim == 3)
      system->time_evolving(_flow_vars.w_var());

    system->time_evolving(_temp_vars.T_var());
    system->time_evolving(_flow_vars.p_var());

    if( _enable_thermo_press_calc )
      system->time_evolving(_p0_var.p0_var());

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesBase<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_species_vars[0])->get_JxW();
    context.get_element_fe(_species_vars[0])->get_phi();
    context.get_element_fe(_species_vars[0])->get_dphi();
    context.get_element_fe(_species_vars[0])->get_xyz();

    context.get_element_fe(_flow_vars.u_var())->get_JxW();
    context.get_element_fe(_flow_vars.u_var())->get_phi();
    context.get_element_fe(_flow_vars.u_var())->get_dphi();
    context.get_element_fe(_flow_vars.u_var())->get_xyz();

    context.get_element_fe(_temp_vars.T_var())->get_JxW();
    context.get_element_fe(_temp_vars.T_var())->get_phi();
    context.get_element_fe(_temp_vars.T_var())->get_dphi();
    context.get_element_fe(_temp_vars.T_var())->get_xyz();

    context.get_element_fe(_flow_vars.p_var())->get_phi();
    context.get_element_fe(_flow_vars.p_var())->get_xyz();

    return;
  }

} // end namespace GRINS
