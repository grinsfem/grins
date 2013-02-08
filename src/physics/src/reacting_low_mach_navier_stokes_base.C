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

// This class
#include "grins/reacting_low_mach_navier_stokes_base.h"

// GRINS
#include "grins/cached_quantities_enum.h"
#include "grins/cea_thermo.h"
#include "grins/cantera_thermo.h"
#include "grins/constant_transport.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"
#include "grins/grins_kinetics.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

namespace GRINS
{
  template<class Mixture>
  ReactingLowMachNavierStokesBase<Mixture>::ReactingLowMachNavierStokesBase(const std::string& physics_name, 
									    const GetPot& input)
    : Physics(physics_name, input),
      _gas_mixture(input),
      _fixed_density( input("Physics/"+reacting_low_mach_navier_stokes+"/fixed_density", false ) ),
      _fixed_rho_value( input("Physics/"+reacting_low_mach_navier_stokes+"/fixed_rho_value", 0.0 ) )
  {
    this->read_input_options(input);
    
    return;
  }

  template<class Mixture>
  ReactingLowMachNavierStokesBase<Mixture>::~ReactingLowMachNavierStokesBase()
  {
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokesBase<Mixture>::read_input_options( const GetPot& input )
  {
    // Read FE family info
    this->_species_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_FE_family", "LAGRANGE") );

    this->_V_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/V_FE_family", "LAGRANGE") );

    this->_P_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/P_FE_family", "LAGRANGE") );

    this->_T_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/T_FE_family", "LAGRANGE") );

    // Read FE family info
    this->_species_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_order", "SECOND") );

    this->_V_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/V_order", "SECOND") );

    this->_P_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/P_order", "FIRST") );

    this->_T_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/T_order", "SECOND") );

    // Read variable naming info
    this->_n_species = input.vector_variable_size("Physics/Chemistry/species");

    _species_var_names.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	/*! \todo Make this prefix string an input option */
	std::string var_name = "w_"+std::string(input( "Physics/Chemistry/species", "DIE!", i ));
	_species_var_names.push_back( var_name );
      }

    this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
    this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
    this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );
    this->_p_var_name = input("Physics/VariableNames/pressure", GRINS::p_var_name_default );
    this->_T_var_name = input("Physics/VariableNames/temperature", GRINS::T_var_name_default );

    // Read thermodynamic state info
    _p0 = input("Physics/"+reacting_low_mach_navier_stokes+"/p0", 0.0 ); /* thermodynamic pressure */

    _enable_thermo_press_calc = input("Physics/"+reacting_low_mach_navier_stokes+"/enable_thermo_press_calc", false );

    if( _enable_thermo_press_calc )
      {
	_p0_var_name = input("Physics/VariableNames/thermo_presure", "p0" );
      }

    // Read gravity vector
    unsigned int g_dim = input.vector_variable_size("Physics/"+reacting_low_mach_navier_stokes+"/g");

    _g(0) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 1 );
  
    if( g_dim == 3)
      _g(2) = input("Physics/"+reacting_low_mach_navier_stokes+"/g", 0.0, 2 );
  
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokesBase<Mixture>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _species_vars.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	_species_vars.push_back( system->add_variable( _species_var_names[i], 
						       this->_species_order, _species_FE_family) );
      }

    _u_var = system->add_variable( _u_var_name, this->_V_order, _V_FE_family);
    _v_var = system->add_variable( _v_var_name, this->_V_order, _V_FE_family);

    if (_dim == 3)
      _w_var = system->add_variable( _w_var_name, this->_V_order, _V_FE_family);
    else
      _w_var = _u_var;

    _p_var = system->add_variable( _p_var_name, this->_P_order, _P_FE_family);
    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);

    /* If we need to compute the thermodynamic pressure, we force this to be a first
       order scalar variable. */
    if( _enable_thermo_press_calc )
      _p0_var = system->add_variable( _p0_var_name, FIRST, SCALAR);

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokesBase<Mixture>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	system->time_evolving( _species_vars[i] );
      }

    system->time_evolving(_u_var);
    system->time_evolving(_v_var);

    if (dim == 3)
      system->time_evolving(_w_var);

    system->time_evolving(_T_var);
    system->time_evolving(_p_var);

    if( _enable_thermo_press_calc )
      system->time_evolving(_p0_var);

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokesBase<Mixture>::init_context( libMesh::FEMContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.element_fe_var[_species_vars[0]]->get_JxW();
    context.element_fe_var[_species_vars[0]]->get_phi();
    context.element_fe_var[_species_vars[0]]->get_dphi();
    context.element_fe_var[_species_vars[0]]->get_xyz();

    context.element_fe_var[_u_var]->get_JxW();
    context.element_fe_var[_u_var]->get_phi();
    context.element_fe_var[_u_var]->get_dphi();
    context.element_fe_var[_u_var]->get_xyz();

    context.element_fe_var[_T_var]->get_JxW();
    context.element_fe_var[_T_var]->get_phi();
    context.element_fe_var[_T_var]->get_dphi();
    context.element_fe_var[_T_var]->get_xyz();

    context.element_fe_var[_p_var]->get_phi();
    context.element_fe_var[_p_var]->get_xyz();

    return;
  }

  //Instantiate
  template class ReactingLowMachNavierStokesBase< IdealGasMixture<CEAThermodynamics,ConstantTransport,Kinetics> >;
#ifdef GRINS_HAVE_CANTERA
  template class ReactingLowMachNavierStokesBase< IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics> >;
  template class ReactingLowMachNavierStokesBase< IdealGasMixture<CanteraThermodynamics,ConstantTransport,CanteraKinetics> >;
  template class ReactingLowMachNavierStokesBase< IdealGasMixture<CEAThermodynamics,ConstantTransport,CanteraKinetics> >;
#endif

} // namespace GRINS
