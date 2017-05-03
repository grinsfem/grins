//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/boussinesq_buoyancy_base.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  BoussinesqBuoyancyBase::BoussinesqBuoyancyBase( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::physics_press_variable_name(input,physics_name))),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _rho(0.0),
      _T_ref(1.0),
      _beta_T(1.0)
  {
    this->read_property(input,
                        "Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/rho_ref",
                        "Density",
                        _rho);

    this->read_property(input,
                        "Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/T_ref",
                        "ReferenceTemperature",
                        _T_ref);


    this->read_property(input,
                        "Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/beta_T",
                        "ThermalExpansionCoeff",
                        _beta_T);

    unsigned int g_dim = input.vector_variable_size("Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/g");

    _g(0) = input("Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/g", 0.0, 1 );

    if( g_dim == 3)
      _g(2) = input("Physics/"+PhysicsNaming::boussinesq_buoyancy()+"/g", 0.0, 2 );

    this->check_var_subdomain_consistency(_flow_vars);
    this->check_var_subdomain_consistency(_press_var);
    this->check_var_subdomain_consistency(_temp_vars);
  }

  void BoussinesqBuoyancyBase::read_property( const GetPot& input,
                                              const std::string& old_option,
                                              const std::string& property,
                                              libMesh::Real& value )
  {
    std::string material = MaterialsParsing::material_name(input,PhysicsNaming::boussinesq_buoyancy());

    // Can't specify both material and rho_ref
    MaterialsParsing::duplicate_input_test(input,
                                           old_option,
                                           "Materials/"+material+"/"+property+"/value" );

    // Deprecated
    if( input.have_variable(old_option) )
      {
        MaterialsParsing::dep_input_warning( old_option,property+"/value" );

        this->set_parameter
          (value, input,
           old_option, 1.0 /*Old default*/);
      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/"+property+"/value" ) )
      {
        this->set_parameter
          (value, input,
           "Materials/"+material+"/"+property+"/value", 0.0 /*default*/);
      }
    // If nothing was set, we default to 1.0. Deprecated, what was I thinking
    else
      {
        this->no_input_warning( input,
                                old_option,
                                material,
                                property );
        this->set_parameter
          (value, input, old_option, 1.0 /*default*/);
      }

    // Make sure density is positive
    if( value <= 0.0 )
      {
        libmesh_error_msg("ERROR: Detected non-positive "+property+"!");
      }
  }

  void BoussinesqBuoyancyBase::no_input_warning( const GetPot& input,
                                                 const std::string& old_option,
                                                 const std::string& material,
                                                 const std::string& property )
  {
    std::string warning;
    if( !MaterialsParsing::have_material(input,PhysicsNaming::boussinesq_buoyancy()) )
      {
        warning = "WARNING: Neither "+old_option+"\n";
        warning += "         nor a material_name was detected in input\n";
        warning +  "         "+property+" is defaulting to 1.0. This behavior is DEPRECATED!\n";
        warning += "         Please update to use Material/MATERIAL_NAME/"+property+"/value\n";
      }
    else
      {
        warning = "WARNING: Neither "+old_option+"\n";
        warning += "         nor Material/"+material+"/"+property+"/value\n";
        warning +  "         "+property+" is defaulting to 1.0. This behavior is DEPRECATED!\n";
        warning += "         Please update to use Material/MATERIAL_NAME/"+property+"/value\n";
      }

    grins_warning(warning);
  }



} // namespace GRINS
