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

// GRINS
#include "chemical_mixture.h"

int test_species( const std::map<std::string,GRINS::Species>& species_name_map,
		  const std::map<GRINS::Species,GRINS::ChemicalSpecies*>& chemical_species,
		  const std::string& species_name,
		  Real molar_mass, Real gas_constant, Real formation_enthalpy, 
		  Real n_tr_dofs, int charge );

int main()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  GRINS::ChemicalMixture chem_mixture( species_str_list );

  const std::map<std::string,GRINS::Species>& species_name_map = chem_mixture.species_name_map();
  const std::map<GRINS::Species,GRINS::ChemicalSpecies*>& chemical_species = chem_mixture.chemical_species();
  const std::vector<GRINS::Species> species_list = chem_mixture.species_list();

  int return_flag = 0;

  for( unsigned int i = 0; i < n_species; i++ )
    {
      if( species_name_map.find( species_str_list[i] )->second != species_list[i] )
	{
	  std::cerr << "Error: species name map and species list ordering mismatch" << std::endl
		    << "species_name_map = " << species_name_map.find( species_str_list[i] )->second
		    << ", species_list = " << species_list[i] << std::endl;
	  return_flag = 1;
	}
    }

  // Check N2 properties
  {
    Real molar_mass = 28.01600;
    return_flag = test_species( species_name_map, chemical_species, "N2",
				molar_mass, GRINS::Constants::R_universal/molar_mass, 0.0, 2.5, 0);
  }
  
  // Check O2 properties
  {
    Real molar_mass = 32.00000;
    return_flag = test_species( species_name_map, chemical_species, "O2",
				molar_mass, GRINS::Constants::R_universal/molar_mass, 0.0, 2.5, 0);
  }

  // Check N properties
  {
    Real molar_mass = 14.00800;
    return_flag = test_species( species_name_map, chemical_species, "N",
				molar_mass, GRINS::Constants::R_universal/molar_mass, 3.3621610000e7, 1.5, 0);
  }

  // Check O properties
  {
    Real molar_mass = 16.00000;
    return_flag = test_species( species_name_map, chemical_species, "O",
				molar_mass, GRINS::Constants::R_universal/molar_mass, 1.5420000000e7, 1.5, 0);
  }

  // Check NO properties
  {
    Real molar_mass = 30.00800;
    return_flag = test_species( species_name_map, chemical_species, "NO",
				molar_mass, GRINS::Constants::R_universal/molar_mass, 2.9961230000e6, 2.5, 0);
  }
  
  return return_flag;
}

int test_species( const std::map<std::string,GRINS::Species>& species_name_map,
		  const std::map<GRINS::Species,GRINS::ChemicalSpecies*>& chemical_species,
		  const std::string& species_name,
		  Real molar_mass, Real gas_constant, Real formation_enthalpy, 
		  Real n_tr_dofs, int charge )
{

  int return_flag = 0;

  GRINS::Species species = species_name_map.find( species_name )->second;
  const GRINS::ChemicalSpecies& chem_species = *(chemical_species.find( species )->second);

  if( chem_species.species() != species_name )
    {
      std::cerr << "Error: Name mismatch for "<< species_name << std::endl
		<< "name = " << chem_species.species() << std::endl;
      return_flag = 1;
    }

  if( chem_species.molar_mass() != molar_mass )
    {
      std::cerr << "Error: Molar mass mismatch for "<< species_name << std::endl
		<< "molar mass = " << chem_species.molar_mass() << std::endl;
      return_flag = 1;
    }

  if( chem_species.gas_constant() != gas_constant )
    {
      std::cerr << "Error: Gas constant mismatch for "<< species_name << std::endl
		<< "gas constant = " << chem_species.gas_constant() << std::endl;
      return_flag = 1;
    }

  if( chem_species.formation_enthalpy() != formation_enthalpy )
    {
      std::cerr << "Error: Formation enthalpy mismatch for "<< species_name << std::endl
		<< "formation enthalpy = " << chem_species.formation_enthalpy() << std::endl;
      return_flag = 1;
    }

  if( chem_species.n_tr_dofs() != n_tr_dofs )
    {
      std::cerr << "Error: Number translational DoFs mismatch for "<< species_name << std::endl
		<< "n_tr_dofs = " << chem_species.n_tr_dofs() << std::endl;
      return_flag = 1;
    }

  if( chem_species.charge() != charge )
    {
      std::cerr << "Error: Charge mismatch for "<< species_name << std::endl
		<< "charge = " << chem_species.charge() << std::endl;
      return_flag = 1;
    }

  return return_flag;
}
