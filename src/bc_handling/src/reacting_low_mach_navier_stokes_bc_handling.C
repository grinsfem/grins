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
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/catalytic_wall.h"

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"

namespace GRINS
{
  ReactingLowMachNavierStokesBCHandling::ReactingLowMachNavierStokesBCHandling( const std::string& physics_name,
										const GetPot& input,
										const ChemicalMixture& chem_mixture )
    : LowMachNavierStokesBCHandling(physics_name,input),
      _n_species( input.vector_variable_size("Physics/Chemistry/species") ),
      _species_var_names(_n_species),
      _species_vars(_n_species),
      _chem_mixture(chem_mixture)
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
	  
	  std::vector<libMesh::Real> species_mass_fracs(n_species_comps);

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
	  std::string reactions_string = "Physics/"+_physics_name+"/wall_catalytic_reactions_"+bc_id_string;
	  if( !input.have_variable(reactions_string) )
	    {
	      std::cerr << "Error: Could not find list of catalytic reactions for boundary id " << bc_id 
			<< std::endl;
	      libmesh_error();
	    }

	  const unsigned int n_reactions = input.vector_variable_size(reactions_string);

	  std::vector<Species> reactants;
	  std::vector<Species> products;

	  // Here, we are assuming 1 reactant and 1 product per reaction
	  reactants.reserve( n_reactions );
	  products.reserve( n_reactions );

	  for( unsigned int r = 0; r < n_reactions; r++ )
	    {
	      std::string reaction = input(reactions_string, "DIE!", r);

	      // First, split each reaction into reactants and products
	      std::vector<std::string> partners;       
	      SplitString(reaction, "->", partners);

	      const std::string& reactant = partners[0];
	      const std::string& product = partners[1];

	      // We currently can only handle reactions of the type R -> P, i.e not R1+R2 -> P, etc.
	      if( partners.size() == 2 )
		{
		  // Parse the reactant and product species and cache
		  const Species& r_species = _chem_mixture.species_name_map().find( reactant )->second;
		  const Species& p_species = _chem_mixture.species_name_map().find( product )->second;

		  // Make sure there's something in there already. Fix for GLIBCXX_DEBUG error.
		  if( _reactant_list.find(bc_id) != _reactant_list.end() )
		    {
		      std::vector<Species>::const_iterator r_it =
			std::find( (_reactant_list.find(bc_id)->second).begin(),
				   (_reactant_list.find(bc_id)->second).end(),
				   r_species );

		      if(  r_it != (_reactant_list.find(bc_id)->second).end() )
			{
			  std::cerr << "Error: Tried adding duplicate reactant " << reactant << " to reactant list."
				    << std::endl;
			  libmesh_error();
			}

		      std::vector<Species>::const_iterator p_it =
			std::find( (_product_list.find(bc_id)->second).begin(),
				   (_product_list.find(bc_id)->second).end(),
				   p_species );

		      if( p_it != (_product_list.find(bc_id)->second).end() )
			{
			  std::cerr << "Error: Tried adding duplicate product " << product << " to product list."
				    << std::endl;
			  libmesh_error();
			}
		    
		    }

		  reactants.push_back( r_species );
		  products.push_back( p_species );

		  // Parse the corresponding catalyticities and cache
		  std::string gamma_r_string = "Physics/"+_physics_name+"/gamma_"+reactant+"_"+bc_id_string;
		  std::string gamma_p_string = "Physics/"+_physics_name+"/gamma_"+product+"_"+bc_id_string;

		  if( !input.have_variable(gamma_r_string) )
		    {
		      std::cout << "Error: Could not find catalyticity for species " << reactant 
				<< ", for boundary " << bc_id << std::endl;
		      libmesh_error();
		    }

		  if( !input.have_variable(gamma_p_string) )
		    {
		      std::cout << "Error: Could not find catalyticity for species " << product 
				<< ", for boundary " << bc_id << std::endl;
		      libmesh_error();
		    }

		  {
		    libMesh::Real gamma_r = input(gamma_r_string, 0.0);

		    if( _catalycities.find(bc_id) == _catalycities.end() )
		      {
			std::map<Species,libMesh::Real> dummy;
			dummy.insert( std::make_pair( r_species, gamma_r ) );
			_catalycities.insert( std::make_pair( bc_id, dummy ) );
		      }
		    else
		      {
			(_catalycities.find(bc_id)->second).insert( std::make_pair( r_species, gamma_r ) );
		      }
		  }

		  {
		    libMesh::Real gamma_p = input(gamma_p_string, 0.0);
			
		    (_catalycities.find(bc_id)->second).insert( std::make_pair( p_species, gamma_p ) );
		  }
		}
	      else
		{
		  std::cerr << "Error: Can currently only handle 1 reactant and 1 product" << std::endl
			    << "in a catalytic reaction." << std::endl
			    << "Found " << partners.size() << " species." << std::endl;
		  libmesh_error();
		}

	    } // end loop over catalytic reactions

	  _reactant_list.insert( std::make_pair( bc_id, reactants ) );
	  _product_list.insert( std::make_pair( bc_id, products ) );
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

  void ReactingLowMachNavierStokesBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    // Call base class
    LowMachNavierStokesBCHandling::init_bc_data(system);

    for( unsigned int s = 0; s < this->_n_species; s++ )
      {
	_species_vars[s] = system.variable_number( _species_var_names[s] );
      }

    // See if we have a catalytic wall
    for( std::map< GRINS::BoundaryID, GRINS::BCType>::const_iterator bc_map = _neumann_bc_map.begin();
	 bc_map != _neumann_bc_map.end(); ++bc_map )
      {
	const BoundaryID bc_id = bc_map->first;
	const BCType bc_type = bc_map->second;

	// Add CatalyticWall for each reactant and product
	if( bc_type == CATALYTIC_WALL )
	  {
	    NBCContainer cont;
	    cont.set_bc_id( bc_id );

	    const std::vector<Species>& reactants = _reactant_list.find(bc_id)->second;
	    const std::vector<Species>& products = _product_list.find(bc_id)->second;

	    /*! \todo  Here we are assuming the same number of reactants and products */
	    libmesh_assert_equal_to( reactants.size(), products.size() );
	    unsigned int n_reactions = reactants.size();
	    
	    for( unsigned int r = 0; r < n_reactions; r++ )
	      {
		/*! \todo  Here we are assuming the same number of reactants and products */
		unsigned int r_species_idx = _chem_mixture.species_list_map().find(reactants[r])->second;
		unsigned int p_species_idx = _chem_mixture.species_list_map().find(products[r])->second;

		libMesh::Real gamma = (_catalycities.find(bc_id)->second).find(reactants[r])->second; 
		
		// -gamma since the reactant is being consumed
		{
		  std::tr1::shared_ptr<NeumannFuncObj> func( new CatalyticWall( _chem_mixture,
										r_species_idx,
										_T_var,
										-gamma ) );
		  
		  VariableIndex var = _species_vars[r_species_idx];
		  
		  cont.add_var_func_pair( var, func );
		}

		// Now products. Using the same gamma as the reactant.
		/*! \todo  Here we are assuming the same number of reactants and products */
		{
		  std::tr1::shared_ptr<NeumannFuncObj> func( new CatalyticWall( _chem_mixture,
										r_species_idx, /* reactant! */
										_T_var,
										gamma ) );
		  
		  VariableIndex var = _species_vars[p_species_idx];
		  
		  cont.add_var_func_pair( var, func );
		}

	      }

	    this->attach_neumann_bound_func( cont );

	  } // end check on CATALYTIC_WALL
      } // end loop over bc_ids

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

	      ConstFunction<libMesh::Number> species_func( this->get_species_bc_value(bc_id,s) );

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
								     const std::vector<libMesh::Real>& species_values )
  {
    _species_bc_values[bc_id] = species_values;
    return;
  }

  libMesh::Real ReactingLowMachNavierStokesBCHandling::get_species_bc_value( GRINS::BoundaryID bc_id, 
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

      case( CATALYTIC_WALL ):
	{
	  libmesh_assert( _reactant_list.find(bc_id) != _reactant_list.end() );
	  libmesh_assert( _product_list.find(bc_id)  != _product_list.end() );

	  const std::vector<Species>& reactants = _reactant_list.find(bc_id)->second;
	  const std::vector<Species>& products = _product_list.find(bc_id)->second;

	  for( std::vector<Species>::const_iterator reactant = reactants.begin();
	       reactant != reactants.end();
	       ++reactant )
	    {
	      /*! \todo We should just cache the map between species and variable number explicitly */
	      std::vector<std::string>::const_iterator it = std::search_n( _species_var_names.begin(),
									   _species_var_names.end(),
									   1,
									   "w_"+_chem_mixture.species_inverse_name_map().find(*reactant)->second );

	      unsigned int species_idx = static_cast<unsigned int>(it - _species_var_names.begin());

	      const VariableIndex var = _species_vars[species_idx];

	      _bound_conds.apply_neumann_normal( context, cache, request_jacobian, var,
						 1.0, 
						 this->get_neumann_bound_func( bc_id, var ) );
	    }

	  for( std::vector<Species>::const_iterator product = products.begin();
	       product != products.end();
	       ++product )
	    {
	      /*! \todo We should just cache the map between species and variable number explicitly */
	      std::vector<std::string>::const_iterator it = std::search_n( _species_var_names.begin(),
									   _species_var_names.end(),
									   1,
									   "w_"+_chem_mixture.species_inverse_name_map().find(*product)->second );

	      unsigned int species_idx = static_cast<unsigned int>(it - _species_var_names.begin());

	      const VariableIndex var = _species_vars[species_idx];

	      _bound_conds.apply_neumann_normal( context, cache, request_jacobian, var,
						 1.0, 
						 this->get_neumann_bound_func( bc_id, var ) );
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
