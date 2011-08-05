//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id:$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "boussinesq_buoyancy.h"

void GRINS::BoussinesqBuoyancy::read_input_options( GetPot& input )
{
  /*
  _rho_ref = ;
  _T_ref = ;
  _beta_T = ;

  _g(1) = ;
  _g(2) = ;
  _g(3) = ;
  */
  return;
}

void GRINS::BoussinesqBuoyancy::init_variables( libMesh::FEMSystem* system )
{
  return;
}

void GRINS::BoussinesqBuoyancy::register_variable_indices(GRINS::VariableMap &global_map)
{
  _u_var = global_map[_u_var_name];
  _v_var = global_map[_v_var_name];
  _w_var = global_map[_w_var_name];
  _T_var = global_map[_T_var_name];

  return;
}

bool GRINS::BoussinesqBuoyancy::element_time_derivative( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
  
  return request_jacobian;
}

void GRINS::BoussinesqBuoyancy::init_context( libMesh::DiffContext &context )
{
  return;
}

bool GRINS::BoussinesqBuoyancy::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext& context,
						      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::element_constraint( bool request_jacobian,
						    libMesh::DiffContext& context,
						    libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::side_constraint( bool request_jacobian,
						 libMesh::DiffContext& context,
						 libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::mass_residual( bool request_jacobian,
					       libMesh::DiffContext& context,
					       libMesh::FEMSystem* system )
{
  return request_jacobian;
}

void GRINS::BoussinesqBuoyancy::build_local_variable_map()
{
  return;
}

