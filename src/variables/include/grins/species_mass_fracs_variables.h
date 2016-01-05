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


#ifndef GRINS_SPECIES_MASS_FRACS_VARIABLES_H
#define GRINS_SPECIES_MASS_FRACS_VARIABLES_H

// C++
#include <vector>

// GRINS
#include "grins/variables_base.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  class SpeciesMassFractionsVariables : public VariablesBase
  {
  public:

    SpeciesMassFractionsVariables( const GetPot& input, const std::string& material_name );
    ~SpeciesMassFractionsVariables(){};

    virtual void init( libMesh::FEMSystem* system );

    unsigned int n_species() const;

    VariableIndex species_var( unsigned int species ) const;

  private:

    SpeciesMassFractionsVariables();

  };

  inline
  unsigned int SpeciesMassFractionsVariables::n_species() const
  {
    // We *must* use the size of _var_names here since that gets populated
    // at construction time.
    return _var_names.size();
  }

  inline
  VariableIndex SpeciesMassFractionsVariables::species_var( unsigned int species ) const
  {
    return _vars[species];
  }

} // end namespace GRINS

#endif // GRINS_SPECIES_MASS_FRACS_VARIABLES_H
