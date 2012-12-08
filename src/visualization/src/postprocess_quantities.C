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

#include "postprocessed_quantities.h"

namespace GRINS
{
  template<class NumericType>
  PostProcessedQuantities<NumericType>::PostProcessedQuantities( const GetPot& input )
    : libMesh::FEMFunctionBase<NumericType>()
  {
    this->build_name_map();

    /* Parse the quantities requested for postprocessing and cache the 
       corresponding enum value */
    unsigned int n_quantities = input.vector_variable_size( "vis-options/output_vars" );
    std::vector<std::string> names(n_quantities);
    _quantities.resize(n_quantities);

    for( unsigned int n = 0; n < n_quantities; n++ )
      {
	names[n] = input("vis-options/output_vars", "DIE!", i);
	
	std::map<std::string, QuantityList>::const_iterator name_it = 
	  _quantity_name_map.find(names[n])->second;
	
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
      }

    return;
  }

  template<class NumericType>
  PostProcessedQuantities<NumericType>::~PostProcessedQuantities()
  {
    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::initialize( const MultiphysicsSystem& system,
							 const libMesh::EquationSystems& equation_systems )
  {
    // Need to cache the MultiphysicsSystem
    _multiphysics_sys = &system;

    libMesh::System& output_system = equation_systems.add_system<libMesh::System>("interior_output");

    // Do sanity check for each of the variables and add variables to the output system as well as
    // cache needed VariableIndex for each of the variables needed from the MultiphysicsSystem
    for( std::vector<QuantityList>::const_iterator it = _quantities.begin();
	 it != _quantities.end(); it++ )
      {
	switch( *it )
	  {
	  case(PERFECT_GAS_DENSITY):
	    {
	      if( !system.has_physics(low_mach_navier_stokes) )
		{
		  std::cerr << "Error: Must have "<< low_mach_navier_stokes 
			    << " enable for perfect gas density calculation."
			    << std::endl;
		  libmesh_error();
		}
	      output_system.add_variable("rho", FIRST);
	    }
	    break;
	    
	  case(MIXTURE_DENSITY):
	    {
	      if( !system.has_physics(reacting_low_mach_navier_stokes) )
		{
		  std::cerr << "Error: Must have "<< reacting_low_mach_navier_stokes 
			    << " enable for mixture gas density calculation."
			    << std::endl;
		  libmesh_error();
		}
	      output_system.add_variable("rho", FIRST);
	    }
	    break;
	    
	  case(PERFECT_GAS_VISCOSITY):
	  case(SPECIES_VISCOSITY):
	  case(MIXTURE_VISCOSITY):
	  case(PERFECT_GAS_THERMAL_CONDUCTIVITY):
	  case(SPECIES_THERMAL_CONDUCTIVITY):
	  case(MIXTURE_THERMAL_CONDUCTIVITY):
	  case(PERFECT_GAS_SPECIFIC_HEAT_P):
	  case(SPECIES_SPECIFIC_HEAT_P):
	  case(MIXTURE_SPECIFIC_HEAT_P):
	  case(PERFECT_GAS_SPECIFIC_HEAT_V):
	  case(SPECIES_SPECIFIC_HEAT_V):
	  case(MIXTURE_SPECIFIC_HEAT_P):
	  case(MOLE_FRACTIONS):
	  case(OMEGA_DOT):
	       break;

	  default:
	    {
	      std::cerr << "Error: Invalid quantity " << *it << std::endl;
	      libmesh_error();
	    }
	  } // end switch

      } // end quantity loop

    return;
  }

  template<class NumericType>
  NumericType PostProcessedQuantities<NumericType>::component( const libMesh::FEMContext& context, 
							       unsigned int i,
							       const libMesh::Point& p,
							       Real /*time*/ )
  {
    // Check if the Elem* is the same between the incoming context and the cached one. If not,
    // reinit the cached MultiphysicsSystem context

    NumericType value = 0.0;

    return value;
  }
    
  template<class NumericType>
  void PostProcessedQuantities<NumericType>::init_context( const libMesh::FEMContext& )
  {
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

    _quantity_name_map["omega_dot"]      = OMEGA_DOT;

    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::operator()( const libMesh::FEMContext&, const libMesh::Point& p,
							 const Real time,
							 libMesh::DenseVector<NumericType>& output )
  {
    for( unsigned int i = 0; i != output.size(); i++ )
      {
	output[i] = this->component(context,i,p,time);
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
  class template PostProcessedQuantities<Number>;

} // namespace GRINS
