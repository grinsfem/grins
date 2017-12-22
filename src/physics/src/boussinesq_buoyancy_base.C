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
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _rho(0.0),
      _T_ref(1.0),
      _beta_T(1.0)
  {
    MaterialsParsing::read_property(input,
                                    "Density",
                                    PhysicsNaming::boussinesq_buoyancy(),
                                    (*this),
                                    _rho);

    MaterialsParsing::read_property(input,
                                    "ReferenceTemperature",
                                    PhysicsNaming::boussinesq_buoyancy(),
                                    (*this),
                                    _T_ref);

    MaterialsParsing::read_property(input,
                                    "ThermalExpansionCoeff",
                                    PhysicsNaming::boussinesq_buoyancy(),
                                    (*this),
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

} // namespace GRINS
