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


// This class
#include "grins/heat_transfer_base.h"

// GRINS
#include "grins_config.h"
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/heat_transfer_macros.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<class K>
  HeatTransferBase<K>::HeatTransferBase( const std::string& physics_name,
                                         const std::string& core_physics_name,
                                         const GetPot& input )
    : Physics(physics_name, input),
      _flow_vars(input,incompressible_navier_stokes),
      _temp_vars(input,heat_transfer),
      _rho(0.0),
      _Cp(1.0),
      _k(input,input("Physics/"+core_physics_name+"/material", "NoMaterial!"))
  {
    this->read_density( core_physics_name, input );

    this->set_parameter
      (this->_Cp, input,
       "Physics/"+core_physics_name+"/Cp", _Cp);

    this->read_input_options(input);

    return;
  }

  template<class K>
  HeatTransferBase<K>::~HeatTransferBase()
  {
    return;
  }

  template<class K>
  void HeatTransferBase<K>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _flow_vars.init(system);
    _temp_vars.init(system);

    return;
  }

  template<class K>
  void HeatTransferBase<K>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_temp_vars.T_var());

    return;
  }

  template<class K>
  void HeatTransferBase<K>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_temp_vars.T_var())->get_JxW();
    context.get_element_fe(_temp_vars.T_var())->get_phi();
    context.get_element_fe(_temp_vars.T_var())->get_dphi();
    context.get_element_fe(_temp_vars.T_var())->get_xyz();

    context.get_side_fe(_temp_vars.T_var())->get_JxW();
    context.get_side_fe(_temp_vars.T_var())->get_phi();
    context.get_side_fe(_temp_vars.T_var())->get_dphi();
    context.get_side_fe(_temp_vars.T_var())->get_xyz();

    return;
  }

  template<class K>
  void HeatTransferBase<K>::register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _k.register_parameter(param_name, param_pointer);
  }

  template<class K>
  void HeatTransferBase<K>::read_density( const std::string& core_physics_name,
                                          const GetPot& input )
  {
    std::string material = input("Physics/"+core_physics_name+"/material", "DIE!");

    // Error if both material/Density and rho are specified
    if( input.have_variable("Physics/"+core_physics_name+"/rho") &&
        input.have_variable("Materials/"+material+"/Density") )
      {
        libmesh_error_msg("ERROR: Can't specify both rho and Density!");
      }

    // It's deprecated to have nothing and default to 1.0
    if( !input.have_variable("Physics/"+core_physics_name+"/rho") &&
        ( !input.have_variable("Physics/"+core_physics_name+"/material") ||
          !input.have_variable("Materials/"+material+"/Density") ) )
      {
        std::string warning = "WARNING: neither Physics/"+core_physics_name+"/rho nor\n";
        warning += "         Physics/"+core_physics_name+"/material options were detected.\n";
        warning += "         We are assuming a density value of 1.0. This is DEPRECATED.\n";
        warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
        grins_warning(warning);

        this->set_parameter
          (this->_rho, input,
           "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
      }

    // It's deprecated to use rho as the density input
    if( input.have_variable("Physics/"+core_physics_name+"/rho") )
      {
        std::string warning = "WARNING: Using input option Physics/"+core_physics_name+"/rho is DEPRECATED.\n";
        warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
        grins_warning(warning);

        this->set_parameter
          (this->_rho, input,
           "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
      }

    // This is the preferred version
    if( input.have_variable("Physics/"+core_physics_name+"/material") &&
        input.have_variable("Materials/"+material+"/Density") )
      {
        this->set_parameter
          (this->_rho, input,
           "Materials/"+material+"/Density/value", 0.0 /*default*/);
      }
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_HEAT_TRANSFER_SUBCLASS(HeatTransferBase);
