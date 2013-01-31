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
#include "grins/reacting_low_mach_navier_stokes.h"

// GRINS
#include "grins/string_utils.h"

namespace GRINS
{
  ReactingLowMachNavierStokesBCHandling::ReactingLowMachNavierStokesBCHandling( const std::string& physics_name,
										const GetPot& input)
    : LowMachNavierStokesBCHandling(physics_name,input),
      _n_species( input.vector_variable_size("Physics/Chemistry/species") ),
      _species_var_names(_n_species),
      _species_vars(_n_species)
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

  void ReactingLowMachNavierStokesBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    // Call base class
    LowMachNavierStokesBCHandling::init_bc_data(system);

    for( unsigned int s = 0; s < this->_n_species; s++ )
      {
	_species_vars[s] = system.variable_number( _species_var_names[s] );
      }

    // See if we have a catalytic wall
    for( std::map< GRINS::BoundaryID, GRINS::BCType>::const_iterator bc_map = _species_bc_map.begin();
	 bc_map != _species_bc_map.end(); ++bc_map )
      {
	const BoundaryID bc_id = bc_map->first;
	const BCType bc_type = bc_map->second;

	// Add CatalyticWall for each reactant and product
	if( bc_type == CATALYTIC_WALL )
	  {
	    libmesh_assert( _catalytic_reactions.find(bc_id) != _catalytic_reactions.end() );
	    const std::vector<std::string>& reactions = _catalytic_reactions.find(bc_id)->second;

	    const unsigned int n_reactions = reactions.size();
	    for( unsigned int r = 0; r < n_reactions; r++ )
	      {
		// First, split each reaction into reactants and products
		std::vector<std::string> partners;       
		SplitString(reactions[r], "->", partners);

		const std::string& reactant = partners[0];
		const std::string& product = partners[1];

		// We currently can only handle reactions of the type R -> P, i.e not R1+R2 -> P, etc.
		if( partners.size() == 2 )
		  {
		    libmesh_assert( _catalycities.find(bc_id) != _catalycities.end() );
		    const libMesh::Real gamma = (_catalycities.find(bc_id)->second)[r];
		    
		    libmesh_not_implemented();
		    //std::tr1::shared_ptr<NeumannFuncObj> func_reactant( new CatalyticWall(chem_mixture,r_index,T_var,gamma) );
		    //std::tr1::shared_ptr<NeumannFuncObj> func_product( new CatalyticWall(chem_mixture,r_index,T_var,-gamma) );

		    // Query to see if there's already a NBCContainer
		    /* If not, instaniate NBCContainer object, populate it,
		       then add it to _neumann_bound_funcs */
		    if( _neumann_bound_funcs.find(bc_id) == _neumann_bound_funcs.end() )
		      {
			NBCContainer container;
			container.set_bc_id( bc_id );
			//container.add_var_func_pair( system.variable_number(reactant), func_reactant );
			//container.add_var_func_pair( system.variable_number(product), func_product );

			_neumann_bound_funcs.insert( std::make_pair( bc_id, container ) );
		      }
		    // If so, add our variable/function pairs to it
		    else
		      {
			NBCContainer& container = _neumann_bound_funcs.find(bc_id)->second;
			
		      }
		  }
		else
		  {
		    std::cerr << "Error: Can currently only handle 1 reactant and 1 product" << std::endl
			      << "in a catalytic reaction." << std::endl
			      << "Found " << partners.size() << " species." << std::endl;
		    libmesh_error();
		  }

	      } // end loop over reaction pairs
	  } // end check on CATALYTIC_WALL
      } // end loop over bc_ids

    return;
  }

  void ReactingLowMachNavierStokesBCHandling::init_bc_types( const GRINS::BoundaryID bc_id, 
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
	  this->set_dirichlet_bc_type( bc_id, bc_type );
	}
	break;

      case(CATALYTIC_WALL):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );

	  // Parse catalytic reactions on this wall

	  // Parse catalycities
	  // Grab gamma value from input
	  std::ostringstream ss;
	  ss << bc_id;
	  libmesh_not_implemented();
	  std::string var_name; // = "Physics/"+_physics_name+"/gamma_"+reactant+"_"+ss.str();
	  
	  if( !input.have_variable( var_name ) )
	    {
	      std::cerr << "Error: Could not find catalyticity "
			<< var_name << std::endl;
	      libmesh_error();
	    }

	}
	break;

      default:
	{
	  LowMachNavierStokesBCHandling::init_bc_types( bc_id, bc_id_string, bc_type, input );
	}
	break;

      } //switch(bc_type)

    return;
  }

  void ReactingLowMachNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
								       libMesh::DofMap& dof_map,
								       GRINS::BoundaryID bc_id,
								       GRINS::BCType bc_type ) const
  {
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
	      std::vector<GRINS::VariableIndex> dbc_vars(1,_species_vars[s]);

	      ConstFunction<Number> species_func( this->get_species_bc_value(bc_id,s) );

	      libMesh::DirichletBoundary species_dbc( dbc_ids, 
						      dbc_vars, 
						      &species_func );
	    
	      dof_map.add_dirichlet_boundary( species_dbc );
	    }
	}
	break;
      
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
								      const GRINS::CachedValues& cache,
								      const bool request_jacobian,
								      const BoundaryID bc_id,
								      const BCType bc_type ) const
  {
    switch( bc_type )
      {
      case( GENERAL_SPECIES ):
      case( CATALYTIC_WALL ):
	{
	  for( std::vector<VariableIndex>::const_iterator var = _species_vars.begin();
	       var != _species_vars.end();
	       ++var )
	    {
	      _bound_conds.apply_neumann_normal( context, cache, request_jacobian, *var,
						 -1.0, 
						 this->get_neumann_bound_func( bc_id, *var ) );
	    }
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
