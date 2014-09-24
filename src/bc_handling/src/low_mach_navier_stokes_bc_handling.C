//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/low_mach_navier_stokes_bc_handling.h"

// GRINS
#include "grins/string_utils.h"

// libMesh
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  LowMachNavierStokesBCHandling::LowMachNavierStokesBCHandling(const std::string& physics_name,
							       const GetPot& input)
    : BCHandlingBase(physics_name),
      _u_var_name( input("Physics/VariableNames/u_velocity", u_var_name_default ) ),
      _v_var_name( input("Physics/VariableNames/v_velocity", v_var_name_default ) ),
      _w_var_name( input("Physics/VariableNames/w_velocity", w_var_name_default ) ),
      _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) )
  {
    std::string id_str = "Physics/"+_physics_name+"/vel_bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/vel_bc_types";
    std::string var_str = "Physics/"+_physics_name+"/vel_bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/vel_bc_values";

    this->read_bc_data( input, id_str, bc_str, var_str, val_str );

    id_str = "Physics/"+_physics_name+"/temp_bc_ids";
    bc_str = "Physics/"+_physics_name+"/temp_bc_types";
    var_str = "Physics/"+_physics_name+"/temp_bc_variables";
    val_str = "Physics/"+_physics_name+"/temp_bc_values";

    this->read_bc_data( input, id_str, bc_str, var_str, val_str );

    return;
  }

  LowMachNavierStokesBCHandling::~LowMachNavierStokesBCHandling()
  {
    return;
  }

  int LowMachNavierStokesBCHandling::string_to_int( const std::string& bc_type ) const
  {
    int bc_type_out;

    if( bc_type == "no_slip" )
      bc_type_out = NO_SLIP;

    else if( bc_type == "zero_x_velocity" )
      bc_type_out = ZERO_X_VELOCITY;

    else if( bc_type == "zero_y_velocity" )
      bc_type_out = ZERO_Y_VELOCITY;

    else if( bc_type == "parabolic_profile" )
      bc_type_out = PARABOLIC_PROFILE;

    else if( bc_type == "prescribed_vel" )
      bc_type_out = PRESCRIBED_VELOCITY;

    else if( bc_type == "general_velocity" )
      bc_type_out = GENERAL_VELOCITY;

    else if( bc_type == "isothermal" )
      bc_type_out = ISOTHERMAL_WALL;

    else if( bc_type == "general_isothermal" )
      bc_type_out = GENERAL_ISOTHERMAL_WALL;
  
    else if( bc_type == "adiabatic" )
      bc_type_out = ADIABATIC_WALL;
  
    else if( bc_type == "prescribed_heat_flux" )
      bc_type_out = PRESCRIBED_HEAT_FLUX;
  
    else if( bc_type == "general_heat_flux" )
      bc_type_out = GENERAL_HEAT_FLUX;

    else
      {
	// Call base class to detect any physics-common boundary conditions
	bc_type_out = BCHandlingBase::string_to_int( bc_type );
      }

    return bc_type_out;
  }

  void LowMachNavierStokesBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _T_var = system.variable_number( _T_var_name );

    return;
  }

  void LowMachNavierStokesBCHandling::init_bc_types( const std::vector<BoundaryID>& bc_ids,
						     const std::string& bc_id_string, 
						     const int bc_type, 
					             const std::string& bc_vars, 
					             const std::string& bc_value, 
						     const GetPot& input )
  {
    switch(bc_type)
      {
      case(NO_SLIP):
	{
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	}
	break;
      case(ZERO_X_VELOCITY):
        {
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	}
	break;
      case(ZERO_Y_VELOCITY):
        {
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	}
	break;
      case(PRESCRIBED_VELOCITY):
	{
          std::vector<std::string> bc_id_strings;
          SplitString (bc_id_string, ":", bc_id_strings);
          libmesh_assert_equal_to(bc_id_strings.size(), bc_ids.size());

          for (int b = 0; b != bc_ids.size(); ++b)
            {
	      this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	
	      /* Force the user to specify 3 velocity components regardless of dimension.
	         This should make it easier to keep things correct if we want to have 
	         2D flow not be in the x-y plane. */
	      int n_vel_comps = input.vector_variable_size("Physics/"+_physics_name+"/bound_vel_"+bc_id_strings[b]);
	      if( n_vel_comps != 3 )
	        {
	          std::cerr << "Error: Must specify 3 velocity components when inputting"
			    << std::endl
			    << "       prescribed velocities. Found " << n_vel_comps
			    << " velocity components."
			    << std::endl;
	          libmesh_error();
	        }
	
	      this->set_dirichlet_bc_value( bc_ids[b], 
					    input("Physics/"+_physics_name+"/bound_vel_"+bc_id_strings[b], 0.0, 0 ),
					    0 );

	      this->set_dirichlet_bc_value( bc_ids[b], 
					    input("Physics/"+_physics_name+"/bound_vel_"+bc_id_strings[b], 0.0, 1 ),
					    1 );

	      this->set_dirichlet_bc_value( bc_ids[b], 
					    input("Physics/"+_physics_name+"/bound_vel_"+bc_id_strings[b], 0.0, 2 ),
					    2 );
	    }
	}
	break;
      case(PARABOLIC_PROFILE):
	{
          std::vector<std::string> bc_id_strings;
          SplitString (bc_id_string, ":", bc_id_strings);
          libmesh_assert_equal_to(bc_id_strings.size(), bc_ids.size());

          for (int b = 0; b != bc_ids.size(); ++b)
            {
	      this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	
	      // Make sure all 6 components are there
	      if( input.vector_variable_size("Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b]) != 6 )
	        {
	          std::cerr << "Error: Must specify 6 components when inputting"
			    << std::endl
			    << "       coefficients for a parabolic profile. Found " 
			    << input.vector_variable_size("Physics/"+_physics_name+"/parabolic_profile_"+bc_id_strings[b])
			    << " components."
			    << std::endl;
	          libmesh_error();
	        }

              libMesh::Real ca = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 0 );
              libMesh::Real cb = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 1 );
              libMesh::Real cc = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 2 );
              libMesh::Real cd = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 3 );
              libMesh::Real ce = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 4 );
              libMesh::Real cf = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_strings[b], 0.0, 5 );

	      std::string var = input( "Physics/"+_physics_name+"/parabolic_profile_var_"+bc_id_strings[b], "DIE!" );
	
	      GRINS::DBCContainer cont;
	      cont.add_var_name( var );
	      cont.add_bc_id( bc_ids[b] );

              std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >
                func( new GRINS::ParabolicProfile(ca,cb,cc,cd,ce,cf) );
	      cont.set_func( func );
	      this->attach_dirichlet_bound_func( cont );
	
	      // Set specified components of Dirichlet data to zero
	      std::string fix_var = input( "Physics/"+_physics_name+"/parabolic_profile_fix_"+bc_id_strings[b], "DIE!" );

	      GRINS::DBCContainer cont_fix;
	      cont_fix.add_var_name( fix_var );
	      cont_fix.add_bc_id( bc_ids[b] );

              std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >
                func_fix( new libMesh::ZeroFunction<libMesh::Number>() );
	      cont_fix.set_func( func_fix );
	      this->attach_dirichlet_bound_func( cont_fix );
	    }
	}
	break;

      case(GENERAL_VELOCITY):
	{
          std::vector<std::string> bc_id_strings;
          SplitString (bc_id_string, ":", bc_id_strings);
          libmesh_assert_equal_to(bc_id_strings.size(), bc_ids.size());

          for (int b = 0; b != bc_ids.size(); ++b)
            {
	      this->set_dirichlet_bc_type( bc_ids[b], bc_type );

              // Set specified components of Dirichlet data to zero.
              // Other component is handled in user BC factory.
	      std::string fix_var = input( "Physics/"+_physics_name+"/general_velocity_fix_"+bc_id_strings[b], "DIE!" );

	      GRINS::DBCContainer cont_fix;
	      cont_fix.add_var_name( fix_var );
	      cont_fix.add_bc_id( bc_ids[b] );

              std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >
                func_fix( new libMesh::ZeroFunction<libMesh::Number>() );
	      cont_fix.set_func( func_fix );
	      this->attach_dirichlet_bound_func( cont_fix );
            }
	}
	break;
      case(ISOTHERMAL_WALL):
	{
          std::vector<std::string> bc_id_strings;
          SplitString (bc_id_string, ":", bc_id_strings);
          libmesh_assert_equal_to(bc_id_strings.size(), bc_ids.size());

          for (int b = 0; b != bc_ids.size(); ++b)
            {
	      this->set_temp_bc_type( bc_ids[b], bc_type );

	      this->set_temp_bc_value( bc_ids[b], input("Physics/"+_physics_name+"/T_wall_"+bc_id_strings[b], 0.0 ) );
            }
	}
	break;
      
      case(GENERAL_ISOTHERMAL_WALL):
	{
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_temp_bc_type( bc_ids[b], bc_type );
	}
	break;

      case(ADIABATIC_WALL):
	{
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_neumann_bc_type( bc_ids[b], bc_type );
	}
	break;
      case(PRESCRIBED_HEAT_FLUX):
	{
          std::vector<std::string> bc_id_strings;
          SplitString (bc_id_string, ":", bc_id_strings);
          libmesh_assert_equal_to(bc_id_strings.size(), bc_ids.size());

          for (int b = 0; b != bc_ids.size(); ++b)
            {
	      this->set_neumann_bc_type( bc_ids[b], bc_type );
	
	      libMesh::RealGradient q_in;
	
	      int num_q_components = input.vector_variable_size("Physics/"+_physics_name+"/q_wall_"+bc_id_strings[b]);
	
	      for( int i = 0; i < num_q_components; i++ )
	        {
	          q_in(i) = input("Physics/"+_physics_name+"/q_wall_"+bc_id_strings[b], 0.0, i );
	        }

	      this->set_neumann_bc_value( bc_ids[b], q_in );
	    }
	}
	break;
      case(GENERAL_HEAT_FLUX):
	{
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_neumann_bc_type( bc_ids[b], bc_type );
	}
	break;

      case(AXISYMMETRIC):
        {
          for (int b = 0; b != bc_ids.size(); ++b)
	    this->set_dirichlet_bc_type( bc_ids[b], bc_type );
	}
	break;

      default:
	{
	  // Call base class to detect any physics-common boundary conditions
	  BCHandlingBase::init_bc_types( bc_ids, bc_id_string, bc_type,
                                         bc_vars, bc_value, input );
	}
      } // End switch(bc_type)
  
    return;
  }

  void LowMachNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
							       libMesh::DofMap& dof_map,
							       BoundaryID bc_id,
							       BCType bc_type ) const
  {
    int dim = system->get_mesh().mesh_dimension();

    VariableIndex T_var = system->variable_number( _T_var_name );
    VariableIndex u_var = system->variable_number( _u_var_name );
    VariableIndex v_var = system->variable_number( _v_var_name );
    VariableIndex w_var = -1;
    if( dim == 3 )
      w_var = system->variable_number( _w_var_name );

    switch( bc_type )
      {
      case(ZERO_X_VELOCITY):
        {
          std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(u_var);

          libMesh::ZeroFunction<libMesh::Number> zero;
	
	  libMesh::DirichletBoundary no_slip_dbc(dbc_ids, 
						 dbc_vars, 
						 &zero );
	
	  dof_map.add_dirichlet_boundary( no_slip_dbc );
        }
        break;
      case(ZERO_Y_VELOCITY):
        {
          std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(v_var);

          libMesh::ZeroFunction<libMesh::Number> zero;
	
	  libMesh::DirichletBoundary no_slip_dbc(dbc_ids, 
						 dbc_vars, 
						 &zero );
	
	  dof_map.add_dirichlet_boundary( no_slip_dbc );
        }
        break;
      case(NO_SLIP):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(u_var);
	  dbc_vars.push_back(v_var);
	  if(dim == 3)
	    dbc_vars.push_back(w_var);
	
          libMesh::ZeroFunction<libMesh::Number> zero;
	
	  libMesh::DirichletBoundary no_slip_dbc(dbc_ids, 
						 dbc_vars, 
						 &zero );
	
	  dof_map.add_dirichlet_boundary( no_slip_dbc );
	}
	break;
      case(PRESCRIBED_VELOCITY):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	
	  // This is inefficient, but it shouldn't matter because
	  // everything gets cached on the libMesh side so it should
	  // only affect performance at startup.
	  {
	    dbc_vars.push_back(u_var);
            libMesh::ConstFunction<libMesh::Number>
              vel_func( this->get_dirichlet_bc_value(bc_id,0) );
	  
	    libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					       dbc_vars, 
					       &vel_func );
	  
	    dof_map.add_dirichlet_boundary( vel_dbc );
	    dbc_vars.clear();
	  }
	
	  {
	    dbc_vars.push_back(v_var);
            libMesh::ConstFunction<libMesh::Number>
              vel_func( this->get_dirichlet_bc_value(bc_id,1) );
	  
	    libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					       dbc_vars, 
					       &vel_func );
	  
	    dof_map.add_dirichlet_boundary( vel_dbc );
	    dbc_vars.clear();
	  }
	  if( dim == 3 )
	    {
	      dbc_vars.push_back(w_var);
              libMesh::ConstFunction<libMesh::Number>
                vel_func( this->get_dirichlet_bc_value(bc_id,2) );
	    
	      libMesh::DirichletBoundary vel_dbc(dbc_ids, 
						 dbc_vars, 
						 &vel_func );

	      dof_map.add_dirichlet_boundary( vel_dbc );
	    }  
	}
	break;
      case(PARABOLIC_PROFILE):
	// This case is handled init_dirichlet_bc_func_objs
	break;

      case(GENERAL_VELOCITY):
	// This case is handled in the BoundaryConditionFactory classes.
	break;

      case(ISOTHERMAL_WALL):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(T_var);
	
          libMesh::ConstFunction<libMesh::Number>
            t_func(this->get_temp_bc_value(bc_id));
	
	  libMesh::DirichletBoundary t_dbc( dbc_ids, dbc_vars, &t_func );
	
	  dof_map.add_dirichlet_boundary( t_dbc );
	}
	break;
      case(GENERAL_ISOTHERMAL_WALL):
	// This case is handled in the BoundaryConditionFactory classes.
	break;

      case(AXISYMMETRIC):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(u_var);
	
          libMesh::ZeroFunction<libMesh::Number> zero;
	
	  libMesh::DirichletBoundary no_slip_dbc( dbc_ids, 
						  dbc_vars, 
						  &zero );
	
	  dof_map.add_dirichlet_boundary( no_slip_dbc );
	}
      break;

      default:
	{
	  std::cerr << "Invalid BCType " << bc_type << std::endl;
	  libmesh_error();
	}
      
      }// end switch

    return;
  }

  void LowMachNavierStokesBCHandling::set_temp_bc_type( BoundaryID bc_id, int bc_type )
  {
    _temp_bc_map.push_back( std::make_pair(bc_id,bc_type) );
    return;
  }

  void LowMachNavierStokesBCHandling::set_temp_bc_value( BoundaryID bc_id,
                                                         libMesh::Real value )
  {
    _T_values[bc_id] = value;
    return;
  }

  libMesh::Real LowMachNavierStokesBCHandling::get_temp_bc_value( BoundaryID bc_id ) const
  {
    return _T_values.find(bc_id)->second;
  }

  void LowMachNavierStokesBCHandling::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
  {
    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::vector<std::pair<BoundaryID,BCType> >::const_iterator it = _dirichlet_bc_map.begin();
         it != _dirichlet_bc_map.end(); it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    for( std::vector<std::pair<BoundaryID,BCType> >::const_iterator it = _temp_bc_map.begin();
        it != _temp_bc_map.end(); it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    return;
  }

} // namespace GRINS
