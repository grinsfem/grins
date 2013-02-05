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
#include "grins/low_mach_navier_stokes_base.h"

// GRINS
#include "grins/constant_viscosity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_conductivity.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  template<class Mu, class SH, class TC>
  LowMachNavierStokesBase<Mu,SH,TC>::LowMachNavierStokesBase(const std::string& physics_name, const GetPot& input)
    : Physics(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  template<class Mu, class SH, class TC>
  LowMachNavierStokesBase<Mu,SH,TC>::~LowMachNavierStokesBase()
  {
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBase<Mu,SH,TC>::read_input_options( const GetPot& input )
  {
    // Read FE info
    this->_V_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/V_FE_family", "LAGRANGE") );

    this->_P_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/P_FE_family", "LAGRANGE") );

    this->_T_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/T_FE_family", "LAGRANGE") );

    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/V_order", "SECOND") );

    this->_P_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/P_order", "FIRST") );

    this->_T_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/T_order", "SECOND") );

    // Read variable naming info
    this->_u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default );
    this->_v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default );
    this->_w_var_name = input("Physics/VariableNames/w_velocity", w_var_name_default );
    this->_p_var_name = input("Physics/VariableNames/pressure", p_var_name_default );
    this->_T_var_name = input("Physics/VariableNames/temperature", T_var_name_default );

    // Read material parameters
    this->_mu.read_input_options( input );
    this->_cp.read_input_options( input );
    this->_k.read_input_options( input );

    // Read thermodynamic state info
    _p0 = input("Physics/"+low_mach_navier_stokes+"/p0", 0.0 ); /* thermodynamic pressure */
    _T0 = input("Physics/"+low_mach_navier_stokes+"/T0", 0.0 ); /* Reference temperature */
    _R  = input("Physics/"+low_mach_navier_stokes+"/R", 0.0 ); /* gas constant */

    if( _R <= 0.0 )
      {
	std::cerr << "=========================================" << std::endl
		  << " Error: Gas constant R must be positive. " << std::endl
		  << " Detected value R = " << _R << std::endl
		  << "=========================================" << std::endl;
	libmesh_error();
      }

    _p0_over_R = _p0/_R;

    _enable_thermo_press_calc = input("Physics/"+low_mach_navier_stokes+"/enable_thermo_press_calc", false );

    if( _enable_thermo_press_calc )
      {
	_p0_var_name = input("Physics/VariableNames/thermo_presure", "p0" );
      }

    // Read gravity vector
    unsigned int g_dim = input.vector_variable_size("Physics/"+low_mach_navier_stokes+"/g");

    _g(0) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 1 );
  
    if( g_dim == 3)
      _g(2) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 2 );
  
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBase<Mu,SH,TC>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

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

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBase<Mu,SH,TC>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

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

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBase<Mu,SH,TC>::init_context( libMesh::FEMContext &context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
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

} // namespace GRINS

// Instantiate
template class GRINS::LowMachNavierStokesBase<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
