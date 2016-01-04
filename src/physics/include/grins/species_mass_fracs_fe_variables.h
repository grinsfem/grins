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
#include "grins/grins_enums.h"
#include "grins/species_mass_fracs_variables.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

namespace GRINS
{

  class SpeciesMassFractionsFEVariables : public SpeciesMassFractionsVariables
  {
  public:

    SpeciesMassFractionsFEVariables( const GetPot& input, const std::string& physics_name );
    ~SpeciesMassFractionsFEVariables(){};

    virtual void init( libMesh::FEMSystem* system );

  protected:

    //! Element type, read from input
    GRINSEnums::FEFamily _species_FE_family;

    //! Element orders, read from input
    GRINSEnums::Order _species_order;

  private:

    SpeciesMassFractionsFEVariables();

  };

} // end namespace GRINS

#endif //GRINS_SPECIES_MASS_FRACS_FE_VARIABLES_H
