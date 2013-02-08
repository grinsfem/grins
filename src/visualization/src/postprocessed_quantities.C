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

#include "grins/postprocessed_quantities.h"

namespace GRINS
{
  template<class NumericType>
  PostProcessedQuantities<NumericType>::PostProcessedQuantities( const GetPot& input )
    : libMesh::FEMFunctionBase<NumericType>(),
      _prev_point(1.0e15,1.0e15,1.0e15) //Initialize to an absurd value
  {
    this->build_name_map();

    /* Parse the quantities requested for postprocessing and cache the 
       corresponding enum value */
    unsigned int n_quantities = input.vector_variable_size( "vis-options/output_vars" );
    std::vector<std::string> names(n_quantities);
    _quantities.resize(n_quantities);

    for( unsigned int n = 0; n < n_quantities; n++ )
      {
	names[n] = input("vis-options/output_vars", "DIE!", n);
	
	typename std::map<std::string, unsigned int>::const_iterator name_it = 
	  _quantity_name_map.find(names[n]);
	
	if( name_it != _quantity_name_map.end() )
	  {
	    _quantities[n] = name_it->second;
	  }
	else
	  {
	    std::cerr << "Error: Invalid name " << names[n] << " for PostProcessedQuantity." 
		      << std::endl;
	    libmesh_error();
	  }

	// Need to cache species names if needed.
	if( names[n] == std::string("mole_fractions") )
	  {
	    unsigned int species_size = input.vector_variable_size( "Physics/Chemistry/species" );
	    _species_names.resize(species_size);
	    
	    for( unsigned int s = 0; s < species_size; s++ )
	      {
		_species_names[s] = input("Physics/Chemistry/species","DIE!",s);
	      }
	  }
      }

    return;
  }

  template<class NumericType>
  PostProcessedQuantities<NumericType>::~PostProcessedQuantities()
  {
    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::initialize( MultiphysicsSystem& system,
							 libMesh::EquationSystems& equation_systems )
  {
    // Only need to initialize if the user requested any output quantities.
    if( !_quantities.empty() )
      {
	// Need to cache the MultiphysicsSystem
	_multiphysics_sys = &system;
 
	libMesh::System& output_system = equation_systems.add_system<libMesh::System>("interior_output");

	// Do sanity check for each of the variables and add variables to the output system as well as
	// cache needed VariableIndex for each of the variables needed from the MultiphysicsSystem
	for( typename std::vector<unsigned int>::const_iterator it = _quantities.begin();
	     it != _quantities.end(); it++ )
	  {
	    this->init_quantities(system,output_system,*it);
	  }

      }
    
    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::init_quantities( const MultiphysicsSystem& multiphysics_system,
							      libMesh::System& output_system,
							      const unsigned int component )
  {
    switch( component )
      {
      case(PERFECT_GAS_DENSITY):
	{
	  if( !multiphysics_system.has_physics(low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< low_mach_navier_stokes 
			<< " enable for perfect gas density calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("rho", FIRST), PERFECT_GAS_DENSITY) );

	  _cache.add_quantity(Cache::PERFECT_GAS_DENSITY);
	}
	break;
	    
      case(MIXTURE_DENSITY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mixture gas density calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("rho", FIRST), MIXTURE_DENSITY) );

	  _cache.add_quantity(Cache::MIXTURE_DENSITY);
	}
	break;
	    
      case(PERFECT_GAS_VISCOSITY):
	{
	  libmesh_not_implemented();
	}
	break;
      case(SPECIES_VISCOSITY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for species viscosity calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("mu_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, SPECIES_VISCOSITY) );
	    }
	  // We need T, p0, and mass fractions too
	  _cache.add_quantity(Cache::TEMPERATURE);
	  _cache.add_quantity(Cache::THERMO_PRESSURE);
	  _cache.add_quantity(Cache::MASS_FRACTIONS);
	  _cache.add_quantity(Cache::SPECIES_VISCOSITY);
	}
	break;

      case(MIXTURE_VISCOSITY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mixture viscosity calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("mu", FIRST), MIXTURE_VISCOSITY) );

	  _cache.add_quantity(Cache::MIXTURE_VISCOSITY);
	}
	break;

      case(PERFECT_GAS_THERMAL_CONDUCTIVITY):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_THERMAL_CONDUCTIVITY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for species thermal conductivity calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("k_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, SPECIES_THERMAL_CONDUCTIVITY) );
	    }

	  _cache.add_quantity(Cache::SPECIES_THERMAL_CONDUCTIVITY);
	}
	break;

      case(MIXTURE_THERMAL_CONDUCTIVITY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mixture thermal conductivity calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("k", FIRST), MIXTURE_THERMAL_CONDUCTIVITY) );

	  _cache.add_quantity(Cache::MIXTURE_THERMAL_CONDUCTIVITY);
	}
	break;

      case(PERFECT_GAS_SPECIFIC_HEAT_P):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_SPECIFIC_HEAT_P):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for species cp calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("cp_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, SPECIES_SPECIFIC_HEAT_P) );
	    }

	  _cache.add_quantity(Cache::SPECIES_SPECIFIC_HEAT_P);
	}
	break;

      case(MIXTURE_SPECIFIC_HEAT_P):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mixture cp calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("cp", FIRST), MIXTURE_SPECIFIC_HEAT_P) );

	  _cache.add_quantity(Cache::MIXTURE_SPECIFIC_HEAT_P);
	}
	break;

      case(PERFECT_GAS_SPECIFIC_HEAT_V):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_SPECIFIC_HEAT_V):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for species cv calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("cv_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, SPECIES_SPECIFIC_HEAT_V) );
	    }

	  _cache.add_quantity(Cache::SPECIES_SPECIFIC_HEAT_V);
	}
	break;

      case(MIXTURE_SPECIFIC_HEAT_V):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mixture cv calculation."
			<< std::endl;
	      libmesh_error();
	    }
	  _quantity_var_map.insert( std::make_pair(output_system.add_variable("cp", FIRST), MIXTURE_SPECIFIC_HEAT_V) );

	  _cache.add_quantity(Cache::MIXTURE_SPECIFIC_HEAT_V);
	}
	break;

      case(MOLE_FRACTIONS):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for mole fraction calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("X_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, MOLE_FRACTIONS) );
	    }

	  _cache.add_quantity(Cache::MOLE_FRACTIONS);
	}
	break;
		
      case(SPECIES_ENTHALPY):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for omega_dot calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("h_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, SPECIES_ENTHALPY) );
	    }

	  // We need T too
	  _cache.add_quantity(Cache::TEMPERATURE);
	  _cache.add_quantity(Cache::SPECIES_ENTHALPY);
	}
	break;

      case(OMEGA_DOT):
	{
	  if( !multiphysics_system.has_physics(reacting_low_mach_navier_stokes) )
	    {
	      std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			<< " enable for omega_dot calculation."
			<< std::endl;
	      libmesh_error();
	    }

	  for( unsigned int s = 0; s < _species_names.size(); s++ )
	    {
	      VariableIndex var = output_system.add_variable("omega_"+_species_names[s], FIRST);
	      _species_var_map.insert( std::make_pair(var, s) );
	      _quantity_var_map.insert( std::make_pair(var, OMEGA_DOT) );
	    }

	  // We need T, p0, and mass fractions too
	  _cache.add_quantity(Cache::TEMPERATURE);
	  _cache.add_quantity(Cache::THERMO_PRESSURE);
	  _cache.add_quantity(Cache::MASS_FRACTIONS);
	  _cache.add_quantity(Cache::MIXTURE_DENSITY);
	  _cache.add_quantity(Cache::MIXTURE_GAS_CONSTANT);
	  _cache.add_quantity(Cache::MOLAR_DENSITIES);
	  _cache.add_quantity(Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY);
	  _cache.add_quantity(Cache::OMEGA_DOT);
	}
	break;

      default:
	{
	  std::cerr << "Error: Invalid quantity " << component << std::endl;
	  libmesh_error();
	}
      } // end switch

    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::update_quantities( libMesh::EquationSystems& equation_systems )
  {
    // Only do the projection if the user actually added any quantities to compute.
    if( !_quantities.empty() )
      {
	libMesh::System& output_system = equation_systems.get_system<libMesh::System>("interior_output");
	output_system.project_solution(this);
      }
    return;
  }
  

  template<class NumericType>
  NumericType PostProcessedQuantities<NumericType>::component( const libMesh::FEMContext& context, 
							       unsigned int component,
							       const libMesh::Point& p,
							       Real /*time*/ )
  {
    // Check if the Elem is the same between the incoming context and the cached one.
    // If not, reinit the cached MultiphysicsSystem context
    if( &(context.get_elem()) != &(_multiphysics_context->get_elem()) )
      {
	_multiphysics_context->pre_fe_reinit(*_multiphysics_sys,&context.get_elem());
	_multiphysics_context->elem_fe_reinit();
      }

    /* Optimization since we expect this function to be called many times with
       the same point. _prev_point initialized to something absurd so this should 
       always be false the first time. */
    if( _prev_point != p )
      {
	_prev_point = p;
	std::vector<libMesh::Point> point_vec(1,p);
	this->_cache.clear();
	_multiphysics_sys->compute_element_cache( *(this->_multiphysics_context), point_vec, this->_cache );
      }

    return this->compute_quantities( component );
  }

  template<class NumericType>
  NumericType PostProcessedQuantities<NumericType>::compute_quantities( const unsigned int component ) const
  {

    NumericType value = 0.0;

    switch( _quantity_var_map.find(component)->second )
      {

      case(PERFECT_GAS_DENSITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::PERFECT_GAS_DENSITY)[0];
	}
	break;
	    
      case(MIXTURE_DENSITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::MIXTURE_DENSITY)[0];
	}
	break;
	    
      case(PERFECT_GAS_VISCOSITY):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_VISCOSITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::SPECIES_VISCOSITY)[0][species];
	}
	break;

      case(MIXTURE_VISCOSITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::MIXTURE_VISCOSITY)[0];
	}
	break;

      case(PERFECT_GAS_THERMAL_CONDUCTIVITY):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_THERMAL_CONDUCTIVITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::SPECIES_THERMAL_CONDUCTIVITY)[0][species];
	}
	break;

      case(MIXTURE_THERMAL_CONDUCTIVITY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY)[0];
	}
	break;

      case(PERFECT_GAS_SPECIFIC_HEAT_P):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_SPECIFIC_HEAT_P):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::SPECIES_SPECIFIC_HEAT_P)[0][species];
	}
	break;

      case(MIXTURE_SPECIFIC_HEAT_P):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[0];
	}
	break;

      case(PERFECT_GAS_SPECIFIC_HEAT_V):
	{
	  libmesh_not_implemented();
	}
	break;

      case(SPECIES_SPECIFIC_HEAT_V):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::SPECIES_SPECIFIC_HEAT_V)[0][species];
	}
	break;

      case(MIXTURE_SPECIFIC_HEAT_V):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  value = this->_cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_V)[0];
	}
	break;

      case(MOLE_FRACTIONS):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::MOLE_FRACTIONS)[0][species];
	}
	break;

      case(SPECIES_ENTHALPY):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::SPECIES_ENTHALPY)[0][species];
	}
	break;

      case(OMEGA_DOT):
	{
	  // Since we only use 1 libMesh::Point, value will always be 0 index of returned vector
	  libmesh_assert( _species_var_map.find(component) != _species_var_map.end() );
	  unsigned int species = _species_var_map.find(component)->second;
	  value = this->_cache.get_cached_vector_values(Cache::OMEGA_DOT)[0][species];
	}
	break;

      default:
	{
	  std::cerr << "Error: Invalid quantity " << _quantity_var_map.find(component)->second << std::endl;
	  libmesh_error();
	}

      } // end switch

    return value;
  }
    
  template<class NumericType>
  void PostProcessedQuantities<NumericType>::init_context( const libMesh::FEMContext& context )
  {
    // Make sure we prepare shape functions for our output variables.
    /*! \todo I believe this is redundant because it's done in the project_vector call. Double check. */
    for( typename std::map<VariableIndex,unsigned int>::const_iterator it = _quantity_var_map.begin();
	 it != _quantity_var_map.end(); it++ )
      {
	libMesh::FEBase* elem_fe = NULL;
	context.get_element_fe( it->first, elem_fe );
	elem_fe->get_phi();
	elem_fe->get_dphi();
	elem_fe->get_JxW();
	elem_fe->get_xyz();

	libMesh::FEBase* side_fe = NULL;
	context.get_side_fe( it->first, side_fe );
	side_fe->get_phi();
	side_fe->get_dphi();
	side_fe->get_JxW();
	side_fe->get_xyz();
      }

    // Create the context we'll be using to compute MultiphysicsSystem quantities
    _multiphysics_context.reset( new libMesh::FEMContext( *_multiphysics_sys ) );
    _multiphysics_sys->init_context(*_multiphysics_context);
    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::build_name_map()
  {
    _quantity_name_map["rho"]            = PERFECT_GAS_DENSITY;
    _quantity_name_map["rho_mix"]        = MIXTURE_DENSITY;

    _quantity_name_map["mu"]             = PERFECT_GAS_VISCOSITY;
    _quantity_name_map["mu_s"]           = SPECIES_VISCOSITY;
    _quantity_name_map["mu_mix"]         = MIXTURE_VISCOSITY;

    _quantity_name_map["k"]              = PERFECT_GAS_THERMAL_CONDUCTIVITY;
    _quantity_name_map["k_s"]            = SPECIES_THERMAL_CONDUCTIVITY;
    _quantity_name_map["k_mix"]          = MIXTURE_THERMAL_CONDUCTIVITY;

    _quantity_name_map["cp"]             = PERFECT_GAS_SPECIFIC_HEAT_P;
    _quantity_name_map["cp_s"]           = SPECIES_SPECIFIC_HEAT_P;
    _quantity_name_map["cp_mix"]         = MIXTURE_SPECIFIC_HEAT_P;

    _quantity_name_map["cv"]             = PERFECT_GAS_SPECIFIC_HEAT_V;
    _quantity_name_map["cv_s"]           = SPECIES_SPECIFIC_HEAT_V;
    _quantity_name_map["cv_mix"]         = MIXTURE_SPECIFIC_HEAT_V;

    _quantity_name_map["mole_fractions"] = MOLE_FRACTIONS;

    _quantity_name_map["h_s"]            = SPECIES_ENTHALPY;
    
    _quantity_name_map["omega_dot"]      = OMEGA_DOT;

    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::operator()( const libMesh::FEMContext& context, const libMesh::Point& p,
							 const Real time,
							 libMesh::DenseVector<NumericType>& output )
  {
    for( unsigned int i = 0; i != output.size(); i++ )
      {
	output(i) = this->component(context,i,p,time);
      }
    return;
  }

  template<class NumericType>
  NumericType PostProcessedQuantities<NumericType>::operator()( const libMesh::FEMContext&, 
								const libMesh::Point&,
								const Real )
  {
    libmesh_error();
    return 0.0; //dummy
  }

  // Instantiate
  template class PostProcessedQuantities<Real>;

} // namespace GRINS
