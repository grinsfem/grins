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


#ifndef GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H
#define GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H

// GRINS
#include "grins/single_fe_type_variable.h"
#include "grins/species_mass_fracs_variables.h"

namespace GRINS
{

  class SpeciesMassFractionsFEVariables : public SingleFETypeVariable,
                                          public SpeciesMassFractionsVariables
  {
  public:

    SpeciesMassFractionsFEVariables( const GetPot& input, const std::string& physics_name,
                                     bool _is_constraint_var = false)
      :  SingleFETypeVariable(input,physics_name,"species_",this->subsection(),"LAGRANGE","SECOND",_is_constraint_var),
         SpeciesMassFractionsVariables(input, MaterialsParsing::material_name(input,physics_name) )
    {}

    ~SpeciesMassFractionsFEVariables(){};

    virtual void init( libMesh::FEMSystem* system )
    { this->default_fe_init(system, _var_names, _vars ); }

  private:

    SpeciesMassFractionsFEVariables();

  };

} // end namespace GRINS

#endif //GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H
