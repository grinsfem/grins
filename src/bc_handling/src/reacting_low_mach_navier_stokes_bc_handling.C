//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "reacting_low_mach_navier_stokes.h"

namespace GRINS
{
  ReactingLowMachNavierStokesBCHandling::ReactingLowMachNavierStokesBCHandling( const std::string& physics_name,
										const GetPot& input)
    : LowMachNavierStokesBCHandling(physics_name,input),
      _n_species( input.vector_variable_size("Physics/Chemistry/species") ),
      _species_var_names(_n_species)
  {

    for( unsigned int s = 0; s < _n_species; s++ )
      {
	/*! \todo Make this prefix string an input option */
	std::string var_name = "w_"+std::string(input( "Physics/Chemistry/species", "DIE!", s ));
	_species_var_names[s] =  var_name;
      }

    std::string id_str = "Physics/"+_physics_name+"/species_bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/species_bc_types";
    
    this->read_bc_data( input, id_str, bc_str );
    
    return;
  }

  ReactingLowMachNavierStokesBCHandling::~ReactingLowMachNavierStokesBCHandling()
  {
    return;
  }

  int ReactingLowMachNavierStokesBCHandling::string_to_int( const std::string& bc_type ) const
  {
    int bc_type_out;

    if( bc_type == "zero_species_flux" )
      {
	bc_type_out = ZERO_SPECIES_FLUX;
      }
    else if( bc_type == "prescribed_species" )
      {
	bc_type_out = PRESCRIBED_SPECIES;
      }
    else if( bc_type == "catalytic_wall" )
      {
	bc_type_out = CATALYTIC_WALL;
      }
    else if( bc_type == "general_species" )
      {
	bc_type_out = GENERAL_SPECIES;
      }
    else
      {
	bc_type_out = LowMachNavierStokesBCHandling::string_to_int( bc_type );
      }

    return bc_type_out;
  }

  void ReactingLowMachNavierStokesBCHandling::init_bc_data( const GRINS::BoundaryID bc_id, 
							    const std::string& bc_id_string, 
							    const int bc_type, 
							    const GetPot& input )
  {
    switch(bc_type)
      {
      case(ZERO_SPECIES_FLUX):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );
	}
	break;
      case(PRESCRIBED_SPECIES):
	{
	  this->set_species_bc_type( bc_id, bc_type );

	  unsigned int n_species_comps = input.vector_variable_size("Physics/"+_physics_name+"/bound_species_"+bc_id_string);

	  if( n_species_comps != _n_species )
	    {
	      std::cerr << "Error: The number of prescribed species values must match" << std::endl
			<< "       the number of species in the simulation." << std::endl
			<< "n_species       = " << _n_species << std::endl
			<< "n_species_comps = " << n_species_comps << std::endl;
	      libmesh_error();
	    }
	  
	  std::vector<Real> species_mass_fracs(n_species_comps);

	  for( unsigned int s = 0; s < n_species_comps; s++ )
	    {
	      species_mass_fracs[s] = input("Physics/"+_physics_name+"/bound_species_"+bc_id_string, -1.0, s );

	      if( (species_mass_fracs[s] > 1.0) ||
		  (species_mass_fracs[s] < 0.0)   )
		{
		  std::cerr << "Error: prescribed species mass fraction must be between 0.0 and 1.0" << std::endl
			    << "w[" << s << "] = " << species_mass_fracs[s] << std::endl;
		  libmesh_error();
		}
	    }

	  this->set_species_bc_values( bc_id, species_mass_fracs );
	}
	break;
      case(GENERAL_SPECIES):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );
	}
	break;
      case(CATALYTIC_WALL):
	{
	  libmesh_not_implemented();
	}
	break;
      default:
	{
	  LowMachNavierStokesBCHandling::init_bc_data( bc_id, bc_id_string, bc_type, input );
	}
      } //switch(bc_type)

    return;
  }

  void ReactingLowMachNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
								       libMesh::DofMap& dof_map,
								       GRINS::BoundaryID bc_id,
								       GRINS::BCType bc_type ) const
  {
    std::vector<VariableIndex> species_vars(_n_species,-1);
    for( unsigned int s = 0; s < _n_species; s++ )
      {
	species_vars[s] = system->variable_number( _species_var_names[s] );
      }

    switch( bc_type )
      {
      case(ZERO_SPECIES_FLUX):
	// Do nothing BC
	break;
      case(PRESCRIBED_SPECIES):
	{
	  std::set<GRINS::BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);

	  for( unsigned int s = 0; s < _n_species; s++ )
	    {
	      std::vector<GRINS::VariableIndex> dbc_vars(1,species_vars[s]);

	      ConstFunction<Number> species_func( this->get_species_bc_value(bc_id,s) );

	      libMesh::DirichletBoundary species_dbc( dbc_ids, 
						      dbc_vars, 
						      &species_func );
	    
	      dof_map.add_dirichlet_boundary( species_dbc );
	    }
	}
	break;
      case(CATALYTIC_WALL):
	{
	  libmesh_not_implemented();
	}
      case(GENERAL_SPECIES):
	// This case is handled in the BoundaryConditionFactory classes.
	break;
      default:
	{
	  LowMachNavierStokesBCHandling::user_init_dirichlet_bcs(system,dof_map,bc_id,bc_type);
	}
      } //switch( bc_type )

    return;
  }

  void ReactingLowMachNavierStokesBCHandling::set_species_bc_type( GRINS::BoundaryID bc_id, int bc_type )
  {
    _species_bc_map[bc_id] = bc_type;
    return;
  }

  void ReactingLowMachNavierStokesBCHandling::set_species_bc_values( GRINS::BoundaryID bc_id, 
								     const std::vector<Real>& species_values )
  {
    _species_bc_values[bc_id] = species_values;
    return;
  }

  Real ReactingLowMachNavierStokesBCHandling::get_species_bc_value( GRINS::BoundaryID bc_id, 
								    unsigned int species ) const
  {
    return (_species_bc_values.find(bc_id)->second)[species];
  }

  void ReactingLowMachNavierStokesBCHandling::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
  {
    LowMachNavierStokesBCHandling::init_dirichlet_bcs(system);

    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::map< GRINS::BoundaryID,GRINS::BCType >::const_iterator it = _species_bc_map.begin();
	 it != _species_bc_map.end();
	 it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    return;
  }

  void ReactingLowMachNavierStokesBCHandling::user_apply_neumann_bcs( libMesh::FEMContext& context,
								      GRINS::VariableIndex var,
								      bool request_jacobian,
								      GRINS::BoundaryID bc_id,
								      GRINS::BCType bc_type ) const
  {
    switch( bc_type )
      {
      case( GENERAL_SPECIES ):
	{
	  _bound_conds.apply_neumann_normal( context, request_jacobian, var, -1.0, 
					     this->get_neumann_bound_func( bc_id, var ) );
	}
	break;
      default:
      {
	std::cerr << "Error: Invalid Neumann BC type for " << _physics_name
		  << std::endl;
	libmesh_error();
      }
      }
    return;
  }

} // namespace GRINS
