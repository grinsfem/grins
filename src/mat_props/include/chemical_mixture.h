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

#ifndef GRINS_CHEMICAL_MIXTURE_H
#define GRINS_CHEMICAL_MIXTURE_H

// C++
#include <vector>
#include <map>
#include <string>

// GRINS
#include "species_enum.h"
#include "chemical_species.h"

namespace GRINS
{
  //! Class storing chemical mixture properties
  /*!
   * This class manages the list of ChemicalSpecies for a requested set
   * of species from input.
   * \todo This should probably be a singleton class, but being lazy for now.
   */
  class ChemicalMixture
  {
  public:
    
    ChemicalMixture( const std::vector<std::string>& species_list );
    ~ChemicalMixture();

    const std::vector<Species>& species_list() const
    { return _species_list; }

    const std::map<Species,ChemicalSpecies*>& chemical_species() const
    { return _chemical_species; }

    const std::map<std::string,Species>& species_name_map() const
    { return _species_name_map; }

  protected:

    void init_species_name_map();

    void read_species_data();
    void read_species_data( std::istream& in );
    
    /*!
     * Skip comment lines in the header of an ASCII
     * text file prefixed with the comment character
     * 'comment_start'.
     * Originally taken from FIN-S.
     */
    void skip_comment_lines( std::istream &in, const char comment_start);

    std::vector<Species> _species_list;
    std::map<Species,ChemicalSpecies*> _chemical_species;
    std::map<std::string,Species> _species_name_map;
    
  private:
    ChemicalMixture();

  };

} //namespace GRINS

#endif //GRINS_CHEMICAL_MIXTURE_H
