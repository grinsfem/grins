//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_HELPER_H
#define GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_HELPER_H

// C++
#include <vector>
#include <string>

// GRINS
#include "grins/var_typedefs.h"
#include <memory>
#include "grins/neumann_bc_abstract.h"
#include "grins/string_utils.h"
#include "grins/fe_variables_base.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  template<typename FunctionType>
  class ParsedFunctionNeumannBCFactoryHelper
  {
  public:

    ParsedFunctionNeumannBCFactoryHelper(){};

    ~ParsedFunctionNeumannBCFactoryHelper(){};

  protected:

    //! Helper function containing common code
    /*! This way, it's easy to add a new flux variable name to be parsed from input. */
    std::shared_ptr<NeumannBCAbstract> build_neumman_func_common( const GetPot& input,
                                                            MultiphysicsSystem& system,
                                                            const FEVariablesBase& fe_var,
                                                            const std::string& flux_input );

    //! Helper function that builds the right BC object depending on the FunctionType
    std::shared_ptr<NeumannBCAbstract> build_parsed_neumann_func(MultiphysicsSystem& system,
                                                           const std::string& expression,
                                                           VariableIndex var_idx );

    //! Helper function that builds the right CompositeBC object depending on the FunctionType
    std::shared_ptr<NeumannBCAbstract>
    build_composite_parsed_neumann_func(MultiphysicsSystem& system,
                                        const std::vector<std::string>& expressions,
                                        const std::vector<VariableIndex>& var_indices );
  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_HELPER_H
