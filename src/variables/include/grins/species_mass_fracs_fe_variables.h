//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H
#define GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H

// GRINS
#include "grins/single_fe_type_variable.h"
#include "grins/variables_parsing.h"

namespace GRINS
{

  class SpeciesMassFractionsFEVariables : public SingleFETypeVariable
  {
  public:

    SpeciesMassFractionsFEVariables( const GetPot& input, const std::string& physics_name,
                                     bool _is_constraint_var = false);

    SpeciesMassFractionsFEVariables( const std::vector<std::string>& var_names,
                                     const std::vector<VariableIndex>& var_indices,
                                     const std::string& prefix,
                                     const std::string& material )
      : SingleFETypeVariable(var_names,var_indices),
        _prefix(prefix),
        _material(material)
    {}

    ~SpeciesMassFractionsFEVariables(){};

    unsigned int n_species() const;

    VariableIndex species( unsigned int species ) const;

    const std::string& material() const
    { return _material; }

    const std::string& prefix() const
    { return _prefix; }

  private:

    SpeciesMassFractionsFEVariables();

    std::string subsection() const
    { return VariablesParsing::species_mass_fractions_section(); }

    std::string _prefix;

    std::string _material;
  };

  inline
  unsigned int SpeciesMassFractionsFEVariables::n_species() const
  {
    // We *must* use the size of _var_names here since that gets populated
    // at construction time.
    return _var_names.size();
  }

  inline
  VariableIndex SpeciesMassFractionsFEVariables::species( unsigned int species ) const
  {
    return _vars[species];
  }

} // end namespace GRINS

#endif //GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H
