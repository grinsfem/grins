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


// This class
#include "grins_config.h"
#include "grins/boussinesq_buoyancy_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  BoussinesqBuoyancyBase::BoussinesqBuoyancyBase( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _T_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") ) ),
      _T_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") ) ),
      _u_var_name( input("Physics/VariableNames/u_velocity", u_var_name_default ) ),
      _v_var_name( input("Physics/VariableNames/v_velocity", v_var_name_default ) ),
      _w_var_name( input("Physics/VariableNames/w_velocity", w_var_name_default ) ),
      _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) ),
      _rho_ref( input("Physics/"+boussinesq_buoyancy+"/rho_ref", 1.0) ),
      _T_ref( input("Physics/"+boussinesq_buoyancy+"/T_ref", 1.0) ),
      _beta_T( input("Physics/"+boussinesq_buoyancy+"/beta_T", 1.0) )
  {
    unsigned int g_dim = input.vector_variable_size("Physics/"+boussinesq_buoyancy+"/g");

    _g(0) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 1 );
  
    if( g_dim == 3)
      _g(2) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 2 );

    return;
  }

  BoussinesqBuoyancyBase::~BoussinesqBuoyancyBase()
  {
    return;
  }

  void BoussinesqBuoyancyBase::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);
 
    // If these are already added, then we just get the index. 
    _u_var = system->add_variable(_u_var_name, _V_order, _V_FE_family );
    _v_var = system->add_variable(_v_var_name, _V_order, _V_FE_family );
    if (_dim == 3)
      _w_var = system->add_variable(_w_var_name, _V_order, _V_FE_family );

    return;
  }

} // namespace GRINS
