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
  PostProcessedQuantities::PostProcessedQuantities( const GetPot& input )
  {
    this->build_name_map();

    /* Parse the quantities requested for postprocessing and cache the 
       corresponding enum value */
    unsigned int n_quantities = input.vector_variable_size( );
    std::vector<std::string> names(n_quantities);
    _quantities.resize(n_quantities);

    for( unsigned int n = 0; n < n_quantities; n++ )
      {
	names[n] = input( , "DIE!", i);
	
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

  PostProcessedQuantities::~PostProcessedQuantities()
  {
    return;
  }

  void PostProcessedQuantities::initialize( const MultiphysicsSystem& system,
					    const libMesh::EquationSystems& equation_systems )
  {
    for( std::vector<QuantityList>::const_iterator it = _quantities.begin();
	 it != _quantities.end(); it++ )
      {
	/* For each quantity, we first do a check that the correct physics is actually present.
	   Then, we add a variable to the output system. */
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
	    }
	    break;
	    
	  } // end switch

      }

    return;
  }

  void PostProcessedQuantities::build_name_map()
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

} // namespace GRINS
