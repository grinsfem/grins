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
#include "grins/solid_mechanics_bc_handling.h"

// libMesh
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  SolidMechanicsBCHandling::SolidMechanicsBCHandling( const std::string& physics_name,
                                                      const GetPot& input )
    : BCHandlingBase(physics_name),
      _disp_vars(input)
  {
    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";
    std::string var_str = "Physics/"+_physics_name+"/bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/bc_values";

    this->read_bc_data( input, id_str, bc_str, var_str, val_str );

    return;
  }

  SolidMechanicsBCHandling::~SolidMechanicsBCHandling()
  {
    return;
  }

  int SolidMechanicsBCHandling::string_to_int( const std::string& bc_type ) const
  {
    int bc_type_out;

    if( bc_type == "pinned" )
      bc_type_out = PINNED;

    else if( bc_type == "constant_displacement" )
      bc_type_out = CONSTANT_DISPLACEMENT;

    else if( bc_type == "roller_x" )
      bc_type_out = ROLLER_X;

    else if( bc_type == "roller_y" )
      bc_type_out = ROLLER_Y;

    else if( bc_type == "roller_z" )
      bc_type_out = ROLLER_Z;

    else if( bc_type == "constant_traction" )
      bc_type_out = CONSTANT_TRACTION;

    else
      {
        // Call base class to detect any physics-common boundary conditions
        bc_type_out = BCHandlingBase::string_to_int( bc_type );
      }

    return bc_type_out;
  }

  void SolidMechanicsBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _disp_vars.init(const_cast<libMesh::FEMSystem*>(&system));

    return;
  }

  void SolidMechanicsBCHandling::init_bc_types( const BoundaryID bc_id,
                                                const std::string& bc_id_string,
                                                const int bc_type,
                                                const std::string& bc_vars,
                                                const std::string& bc_value,
                                                const GetPot& input )
  {
    switch(bc_type)
      {
      case(PINNED):
        {
          this->set_dirichlet_bc_type( bc_id, bc_type );
        }
        break;

      case(CONSTANT_DISPLACEMENT):
        {
          this->set_dirichlet_bc_type( bc_id, bc_type );

          int n_disp_comps = input.vector_variable_size("Physics/"+_physics_name+"/displacement_"+bc_id_string);

          if( _disp_vars.have_v() )
            {
              if( n_disp_comps < 2 )
                {
                  std::cerr << "Error: Must specify at least 2 displacement components for 2-D problem." << std::endl;
                  libmesh_error();
                }
            }

          if( _disp_vars.have_w() )
            {
              if( n_disp_comps < 3 )
                {
                  std::cerr << "Error: Must specify 3 displacement components for 3-D problem." << std::endl;
                  libmesh_error();
                }
            }

          for( int i = 0; i < n_disp_comps; i++ )
            {
              this->set_dirichlet_bc_value( bc_id,
                                            input("Physics/"+_physics_name+"/displacement_"+bc_id_string, 0.0, i ),
                                            i );
            }
        }
        break;

      case(ROLLER_X):
      case(ROLLER_Y):
      case(ROLLER_Z):
        {
          this->set_dirichlet_bc_type( bc_id, bc_type );
        }
        break;

      case(CONSTANT_TRACTION):
        {
          this->set_neumann_bc_type( bc_id, bc_type );

          libMesh::Gradient t_in;

          int num_t_components = input.vector_variable_size("Physics/"+_physics_name+"/traction_"+bc_id_string);

          if( num_t_components < 3 )
            {
              std::cerr << "Error: Expecting 3 components when specifying traction!" << std::endl
                        << "       Found " << num_t_components << " components." << std::endl;
              libmesh_error();
            }

          for( int i = 0; i < num_t_components; i++ )
            {
              t_in(i) = input("Physics/"+_physics_name+"/traction_"+bc_id_string, 0.0, i );
            }

          this->set_neumann_bc_value( bc_id, t_in );
        }
        break;

      default:
        {
          // Call base class to detect any physics-common boundary conditions
          BCHandlingBase::init_bc_types( bc_id, bc_id_string, bc_type,
                                         bc_vars, bc_value, input );
        }

      }// End switch(bc_type)

    return;
  }

  void SolidMechanicsBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
                                                          libMesh::DofMap& dof_map,
                                                          BoundaryID bc_id,
                                                          BCType bc_type ) const
  {
    VariableIndex u_var = _disp_vars.u_var();
    libmesh_assert(system->has_variable(_disp_vars.u_var_name()));

    VariableIndex v_var;
    if( _disp_vars.have_v() )
      {
        v_var = _disp_vars.v_var();
        libmesh_assert(system->has_variable(_disp_vars.v_var_name()));
      }

    VariableIndex w_var;
    if( _disp_vars.have_w() )
      {
        w_var = _disp_vars.w_var();
        libmesh_assert(system->has_variable(_disp_vars.w_var_name()));
      }

    switch( bc_type )
      {
      case(PINNED):
        {
          std::set<BoundaryID> dbc_ids;
          dbc_ids.insert(bc_id);

          std::vector<VariableIndex> dbc_vars;
          dbc_vars.push_back(u_var);

          if( _disp_vars.have_v() )
            dbc_vars.push_back(v_var);

          if( _disp_vars.have_w() )
            dbc_vars.push_back(w_var);

          libMesh::ZeroFunction<libMesh::Number> zero;

          libMesh::DirichletBoundary no_slip_dbc(dbc_ids,
                                                 dbc_vars,
                                                 &zero );

          dof_map.add_dirichlet_boundary( no_slip_dbc );
        }
        break;

      case(CONSTANT_DISPLACEMENT):
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
              disp_func( this->get_dirichlet_bc_value(bc_id,0) );

            libMesh::DirichletBoundary disp_dbc(dbc_ids,
                                                dbc_vars,
                                                &disp_func );

            dof_map.add_dirichlet_boundary( disp_dbc );
            dbc_vars.clear();
          }

          if( _disp_vars.have_v() )
            {
              dbc_vars.push_back(v_var);
              libMesh::ConstFunction<libMesh::Number>
                disp_func( this->get_dirichlet_bc_value(bc_id,1) );

              libMesh::DirichletBoundary disp_dbc(dbc_ids,
                                                  dbc_vars,
                                                  &disp_func );

              dof_map.add_dirichlet_boundary( disp_dbc );
              dbc_vars.clear();
            }

          if( _disp_vars.have_w() )
            {
              dbc_vars.push_back(w_var);
              libMesh::ConstFunction<libMesh::Number>
                disp_func( this->get_dirichlet_bc_value(bc_id,2) );

              libMesh::DirichletBoundary disp_dbc(dbc_ids,
                                                  dbc_vars,
                                                  &disp_func );

              dof_map.add_dirichlet_boundary( disp_dbc );
            }
        }
        break;

        // Roller is free to move in the x-direction, so pin y and z-directions
      case(ROLLER_X):
        {
          std::set<BoundaryID> dbc_ids;
          dbc_ids.insert(bc_id);

          std::vector<VariableIndex> dbc_vars;

          if( _disp_vars.have_v() )
            dbc_vars.push_back(v_var);

          if( _disp_vars.have_w() )
            dbc_vars.push_back(w_var);

          libMesh::ZeroFunction<libMesh::Number> zero;

          libMesh::DirichletBoundary no_slip_dbc(dbc_ids,
                                                 dbc_vars,
                                                 &zero );

          dof_map.add_dirichlet_boundary( no_slip_dbc );
        }
        break;

        // Roller is free to move in the y-direction, so pin x and z-directions
      case(ROLLER_Y):
        {
          std::set<BoundaryID> dbc_ids;
          dbc_ids.insert(bc_id);

          std::vector<VariableIndex> dbc_vars;
          dbc_vars.push_back(u_var);

          if( _disp_vars.have_w() )
            dbc_vars.push_back(w_var);

          libMesh::ZeroFunction<libMesh::Number> zero;

          libMesh::DirichletBoundary no_slip_dbc(dbc_ids,
                                                 dbc_vars,
                                                 &zero );

          dof_map.add_dirichlet_boundary( no_slip_dbc );
        }
        break;

        // Roller is free to move in the z-direction, so pin x and y-directions
      case(ROLLER_Z):
        {
          std::set<BoundaryID> dbc_ids;
          dbc_ids.insert(bc_id);

          std::vector<VariableIndex> dbc_vars;
          dbc_vars.push_back(u_var);
          dbc_vars.push_back(v_var);

          libMesh::ZeroFunction<libMesh::Number> zero;

          libMesh::DirichletBoundary no_slip_dbc(dbc_ids,
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

  void SolidMechanicsBCHandling::user_apply_neumann_bcs( AssemblyContext& context,
                                                         const GRINS::CachedValues& /*cache*/,
                                                         const bool /*request_jacobian*/,
                                                         const BoundaryID bc_id,
                                                         const BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(CONSTANT_TRACTION):
        {
          const libMesh::Point& traction = this->get_neumann_bc_value(bc_id);

          _bound_conds.apply_neumann_normal( context, _disp_vars.u_var(), 1.0, traction(0) );

          if( _disp_vars.have_v() )
            _bound_conds.apply_neumann_normal( context, _disp_vars.v_var(), 1.0, traction(1) );

          if( _disp_vars.have_w() )
            _bound_conds.apply_neumann_normal( context, _disp_vars.w_var(), 1.0, traction(2) );
        }
        break;

      default:
        {
          std::cerr << "Error: Invalid Neumann BC type for " << _physics_name
                    << std::endl;
          libmesh_error();
        }

      } // switch( bc_type )

    return;
  }

} // end namespace GRINS
