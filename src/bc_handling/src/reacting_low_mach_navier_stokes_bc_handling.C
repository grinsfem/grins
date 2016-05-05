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
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/gas_recombination_catalytic_wall.h"
#include "grins/gas_solid_catalytic_wall.h"
#include "grins/constant_catalycity.h"
#include "grins/arrhenius_catalycity.h"
#include "grins/power_law_catalycity.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"
#include "grins/catalycity_factory_old_style_base.h"

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"

namespace GRINS
{
  template<typename Chemistry>
  ReactingLowMachNavierStokesBCHandling<Chemistry>::ReactingLowMachNavierStokesBCHandling( const std::string& physics_name,
                                                                                           const GetPot& input,
                                                                                           const Chemistry& chemistry )
    : LowMachNavierStokesBCHandling(physics_name,input),
      _species_vars(input,MaterialsParsing::material_name(input,PhysicsNaming::reacting_low_mach_navier_stokes())),
      _n_species(_species_vars.n_species()),
      _chemistry(chemistry)
  {
    std::string id_str = "Physics/"+_physics_name+"/species_bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/species_bc_types";
    std::string var_str = "Physics/"+_physics_name+"/species_bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/species_bc_values";
    
    this->read_bc_data( input, id_str, bc_str, var_str, val_str );
    
    return;
  }

  template<typename Chemistry>
  ReactingLowMachNavierStokesBCHandling<Chemistry>::~ReactingLowMachNavierStokesBCHandling()
  {
    return;
  }

  template<typename Chemistry>
  int ReactingLowMachNavierStokesBCHandling<Chemistry>::string_to_int( const std::string& bc_type ) const
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
    else if( bc_type == "prescribed_mole_fracs" )
      {
	bc_type_out = PRESCRIBED_MOLE_FRACTIONS;
      }
    else if( bc_type == "gas_recombination_catalytic_wall" )
      {
        bc_type_out = GAS_RECOMBINATION_CATALYTIC_WALL;
      }
    else if( bc_type == "gas_solid_catalytic_wall" )
      {
        bc_type_out = GAS_SOLID_CATALYTIC_WALL;
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

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::init_bc_types( const GRINS::BoundaryID bc_id, 
							     const std::string& bc_id_string, 
							     const int bc_type, 
					                     const std::string& bc_vars, 
							     const std::string& bc_value, 
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

      case(PRESCRIBED_MOLE_FRACTIONS):
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
	  
	  std::vector<libMesh::Real> species_mole_fracs(n_species_comps);

	  for( unsigned int s = 0; s < n_species_comps; s++ )
	    {
	      species_mole_fracs[s] = input("Physics/"+_physics_name+"/bound_species_"+bc_id_string, -1.0, s );

	      if( (species_mole_fracs[s] > 1.0) ||
		  (species_mole_fracs[s] < 0.0)   )
		{
		  std::cerr << "Error: prescribed species mole fraction must be between 0.0 and 1.0" << std::endl
			    << "w[" << s << "] = " << species_mole_fracs[s] << std::endl;
		  libmesh_error();
		}
	    }

          // Compute M
          libMesh::Real M = 0.0;
          for( unsigned int s = 0; s < n_species_comps; s++ )
	    {
              M += species_mole_fracs[s]*_chemistry.M(s);
            }

          // Convert mole fractions to mass fractions
          std::vector<libMesh::Real> species_mass_fracs(n_species_comps);
          for( unsigned int s = 0; s < n_species_comps; s++ )
	    {
              species_mass_fracs[s] = species_mole_fracs[s]*_chemistry.M(s)/M;
            }

	  this->set_species_bc_values( bc_id, species_mass_fracs );
	}
	break;

      case(GENERAL_SPECIES):
	{
	  this->set_dirichlet_bc_type( bc_id, bc_type );
	}
	break;

      case(GAS_RECOMBINATION_CATALYTIC_WALL):
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

          for( unsigned int r = 0; r < n_reactions; r++ )
	    {
              std::string reaction = input(reactions_string, "DIE!", r);

	      // First, split each reaction into reactants and products
	      std::vector<std::string> partners;
	      StringUtilities::split_string(reaction, "->", partners);

	      const std::string& reactant = partners[0];
	      const std::string& product = partners[1];

	      /*! \todo We currently can only handle reactions of the type R -> P, i.e not R1+R2 -> P, etc. */
	      if( partners.size() == 2 )
		{
		  /* ------------- Grab the reactant and product species indices ------------- */
		  const unsigned int r_species = _chemistry.species_index( reactant );
		  const unsigned int p_species = _chemistry.species_index( product );

                  /* ------------- Parse and construct the corresponding catalyticities ------------- */

                  // These are temporary and will be cloned, so let them be destroyed when we're done
                  libMesh::UniquePtr<CatalycityBase> gamma_r =
                    this->build_catalycities( input, reactant, bc_id_string );

                  /* ------------- Now cache the CatalyticWall functions to init later ------------- */
                  libmesh_assert( gamma_r );

                  SharedPtr<CatalyticWallBase<Chemistry> > wall_ptr( new GasRecombinationCatalyticWall<Chemistry>( _chemistry, *gamma_r, r_species, p_species ) );

                  _catalytic_walls.insert( std::make_pair(bc_id, wall_ptr ) );

                } // if( partners.size() == 2 )
              else
                {
                  std::cerr << "Error: Can currently only handle 1 reactant and 1 product" << std::endl
                            << "in a catalytic reaction." << std::endl
                            << "Found " << partners.size() << " species." << std::endl;
                  libmesh_error();
                }

            } // loop over reactions
        }
        break;

        case(GAS_SOLID_CATALYTIC_WALL):
          {
            this->set_neumann_bc_type( bc_id, bc_type );

            // Parse catalytic reactions on this wall
            std::string reactions_string = "Physics/"+_physics_name+"/wall_gas_solid_reactions_"+bc_id_string;
            if( !input.have_variable(reactions_string) )
              {
                std::cerr << "Error: Could not find list of gas-solid catalytic reactions for boundary id "
                          << bc_id
                          << std::endl;
                libmesh_error();
              }

            const unsigned int n_reactions = input.vector_variable_size(reactions_string);

            for( unsigned int r = 0; r < n_reactions; r++ )
              {
                std::string reaction = input(reactions_string, "DIE!", r);

                /* We are expecting reactions of the form
                   X+Y(s)->Z  or
                   Y(s)+X->X
                So, first we'll split on the "->", then split the reactants up and
                figure out which is the gas species and which is the solid species. */

                std::vector<std::string> partners;
                StringUtilities::split_string(reaction, "->", partners);

                const std::string pre_split_reactants = partners[0];
                const std::string& product = partners[1];

                std::vector<std::string> split_reactants;
                StringUtilities::split_string(pre_split_reactants, "+", split_reactants);

                // We can only handle two reactants currently
                if( split_reactants.size() != 2 )
                  {
                    std::cerr << "Error: Currently, GasSolidCatalyticWall boundary condition only supports"
                              << std::endl
                              << "       reactions of the form X+Y(s)->Z or Y(s)+X->X. Found "
                              << split_reactants.size() << " reactants." << std::endl;
                    libmesh_error();
                  }

                std::string gas_reactant;
                std::string solid_reactant;
                // Check if the first reactant is the solid one
                if( split_reactants[0].find("(s)") == split_reactants[0].npos )
                  {
                    // If not found, check the second entry
                    if( split_reactants[1].find("(s)") == split_reactants[1].npos )
                      {
                        std::cerr << "Error: could not find solid reactant for GasSolidCatalyticWall" << std::endl
                                  << "       boundary condition. Found reactants " << split_reactants[0]
                                  << ", " << split_reactants[1] << std::endl;
                        libmesh_error();
                      }
                    else
                      {
                        gas_reactant = split_reactants[0];
                        solid_reactant = split_reactants[1].substr(0,split_reactants[1].find("(s)"));
                      }
                  }
                // Found (s) in the first entry
                else
                  {
                    // Check that there's not 2 solid reactants
                    if( split_reactants[1].find("(s)") != split_reactants[1].npos )
                      {
                        std::cerr << "Error: can have only one solid reactant for GasSolidCatalyticWall" << std::endl
                                  << "       boundary condition. Found reactants " << split_reactants[0]
                                  << ", " << split_reactants[1] << std::endl;
                        libmesh_error();
                      }

                    gas_reactant = split_reactants[1];
                    solid_reactant = split_reactants[0].substr(0,split_reactants[0].find("(s)"));
                  }

                /* Now we have the gas reactant, the solid reactant, and the gas product strings.
                   Next we grab the species indices, build the catalycity, then build the
                   CatalyticWallBase object. */
                const unsigned int rg_species = _chemistry.species_index( gas_reactant );
                const unsigned int rs_species = _chemistry.species_index( solid_reactant );
                const unsigned int p_species  = _chemistry.species_index( product );

                // This is temporary and will be cloned, so let it be destroyed when we're done
                libMesh::UniquePtr<CatalycityBase> gamma_r =
                  this->build_catalycities( input, gas_reactant, bc_id_string );

                libmesh_assert( gamma_r );

                SharedPtr<CatalyticWallBase<Chemistry> > wall_ptr( new GasSolidCatalyticWall<Chemistry>( _chemistry, *gamma_r, rg_species, rs_species, p_species ) );

                _catalytic_walls.insert( std::make_pair(bc_id, wall_ptr ) );

              } // loop over reactions
          }
        break;

      default:
	{
	  LowMachNavierStokesBCHandling::init_bc_types( bc_id, bc_id_string, bc_type,
                                                        bc_vars, bc_value, input );
	}
	break;

      } //switch(bc_type)

    return;
  }

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::init_bc_data( const libMesh::FEMSystem& system )
  {
    // Call base class
    LowMachNavierStokesBCHandling::init_bc_data(system);

    _species_vars.init_vars( const_cast<libMesh::FEMSystem*>(&system) );

    // See if we have a catalytic wall and initialize them if we do
    for( std::map< GRINS::BoundaryID, GRINS::BCType>::const_iterator bc_map = _neumann_bc_map.begin();
	 bc_map != _neumann_bc_map.end(); ++bc_map )
      {
	const BoundaryID bc_id = bc_map->first;
	const BCType bc_type = bc_map->second;

        if( bc_type == GAS_RECOMBINATION_CATALYTIC_WALL ||
            bc_type == GAS_SOLID_CATALYTIC_WALL )
          {
            typedef typename std::multimap<BoundaryID, SharedPtr<CatalyticWallBase<Chemistry> > >::iterator it_type;

            std::pair< it_type, it_type > it_range = _catalytic_walls.equal_range( bc_id );

            for( it_type it = it_range.first; it != it_range.second; ++it )
              {
                (it->second)->init(system);

                if( this->is_axisymmetric() )
                  {
                    (it->second)->set_axisymmetric(true);
                  }
              }
          } // if( bc_type == GAS_RECOMBINATION_CATALYTIC_WALL )

      } // end loop over bc_ids

    return;
  }

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
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
      case(PRESCRIBED_MOLE_FRACTIONS):
	{
	  std::set<GRINS::BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);

	  for( unsigned int s = 0; s < _n_species; s++ )
	    {
	      std::vector<GRINS::VariableIndex> dbc_vars(1,_species_vars.species(s));

              libMesh::ConstFunction<libMesh::Number> species_func( this->get_species_bc_value(bc_id,s) );

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

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::set_species_bc_type( GRINS::BoundaryID bc_id, int bc_type )
  {
    _species_bc_map.push_back( std::make_pair(bc_id,bc_type) );
    return;
  }

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::set_species_bc_values( GRINS::BoundaryID bc_id, 
                                                                                const std::vector<libMesh::Real>& species_values )
  {
    _species_bc_values[bc_id] = species_values;
    return;
  }

  template<typename Chemistry>
  libMesh::Real ReactingLowMachNavierStokesBCHandling<Chemistry>::get_species_bc_value( GRINS::BoundaryID bc_id, 
									     unsigned int species ) const
  {
    return (_species_bc_values.find(bc_id)->second)[species];
  }

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
  {
    LowMachNavierStokesBCHandling::init_dirichlet_bcs(system);

    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::vector<std::pair<BoundaryID,BCType> >::const_iterator it = _species_bc_map.begin();
	 it != _species_bc_map.end();
	 it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    return;
  }

  template<typename Chemistry>
  void ReactingLowMachNavierStokesBCHandling<Chemistry>::user_apply_neumann_bcs( AssemblyContext& context,
								      const GRINS::CachedValues& cache,
								      const bool request_jacobian,
								      const BoundaryID bc_id,
								      const BCType bc_type ) const
  {
    switch( bc_type )
      {
      case( GENERAL_SPECIES ):
	{
	  for( unsigned int s = 0; s < this->_n_species; s++ )
	    {
              unsigned int var = _species_vars.species(s);

              if( this->is_axisymmetric() )
                {
                  _bound_conds.apply_neumann_normal_axisymmetric( context, cache,
                                                                  request_jacobian, var, -1.0,
                                                                  this->get_neumann_bound_func( bc_id, var ) );
                }
              else
                {
                  _bound_conds.apply_neumann_normal( context, cache,
                                                     request_jacobian, var, 1.0,
                                                     this->get_neumann_bound_func( bc_id, var ) );
                }
	    }
	}
	break;

      case( GAS_RECOMBINATION_CATALYTIC_WALL ):
      case( GAS_SOLID_CATALYTIC_WALL ):
        {
          typedef typename std::multimap<BoundaryID, SharedPtr<CatalyticWallBase<Chemistry> > >::const_iterator it_type;

          std::pair< it_type, it_type > it_range = _catalytic_walls.equal_range( bc_id );

          for( it_type it = it_range.first; it != it_range.second; ++it )
              {
                (it->second)->apply_fluxes(context, cache, request_jacobian );
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

  template<typename Chemistry>
  libMesh::UniquePtr<CatalycityBase>
  ReactingLowMachNavierStokesBCHandling<Chemistry>::build_catalycities( const GetPot& input,
                                                                        const std::string& reactant,
                                                                        const std::string& bc_id_string )
  {
    CatalycityFactoryOldStyleBase::set_getpot(input);
    CatalycityFactoryOldStyleBase::set_section("Physics/"+_physics_name);
    CatalycityFactoryOldStyleBase::set_reactant(reactant);
    CatalycityFactoryOldStyleBase::set_bc_id(bc_id_string);

    std::string catalycity_type = input("Physics/"+_physics_name+"/gamma_"+reactant+"_"+bc_id_string+"_type", "none");

    // "_old_style" added at the end since this parsing style is the "old style".
    return CatalycityFactoryOldStyleBase::build(catalycity_type+"_old_style");
  }

  template<typename Chemistry>
  CatalyticWallBase<Chemistry>* ReactingLowMachNavierStokesBCHandling<Chemistry>::get_catalytic_wall( const BoundaryID bc_id )
  {
    return (_catalytic_walls.find(bc_id)->second).get();
  }

} // namespace GRINS
