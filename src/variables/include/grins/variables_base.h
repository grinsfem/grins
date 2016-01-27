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


#ifndef GRINS_VARIABLES_BASE_H
#define GRINS_VARIABLES_BASE_H

// C++
#include <vector>
#include <string>

// GRINS
#include "grins/var_typedefs.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  class VariablesBase
  {
  public:

    VariablesBase(){};

    ~VariablesBase(){};

  protected:

    //! Default method for init'ing variables
    /*! This method assumes that the variable has already been
        added to the System elsewhere and we're just grabbing the
        variable number from the system. We attempt to grab
        a variable number for every entry in _var_names. */
    void default_var_init( libMesh::FEMSystem* system );

    std::vector<VariableIndex> _vars;

    std::vector<std::string> _var_names;

  };

} // end namespace GRINS

#endif // GRINS_VARIABLES_BASE_H
