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

#ifndef GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_H
#define GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_H

// GRINS
#include "grins/neumann_bc_factory_abstract.h"
#include "grins/parsed_function_neumann_bc_factory_helper.h"
#include "grins/boundary_condition_names.h"

namespace GRINS
{
  template<typename FunctionType>
  class ParsedFunctionNeumannBCFactory : public NeumannBCFactoryAbstract,
                                         public ParsedFunctionNeumannBCFactoryHelper<FunctionType>
  {
  public:

    ParsedFunctionNeumannBCFactory( const std::string& bc_type_name )
      : NeumannBCFactoryAbstract(bc_type_name),
        ParsedFunctionNeumannBCFactoryHelper<FunctionType>()
    {}

    ~ParsedFunctionNeumannBCFactory(){};

  protected:

    virtual std::shared_ptr<NeumannBCAbstract>
    build_neumann_func( const GetPot& input,
                        MultiphysicsSystem& system,
                        const FEVariablesBase& fe_var,
                        const std::string& section )
    {
      std::string flux_input = this->flux_input(section);

      // Make sure flux input specified and consistent with var_names size
      this->check_for_flux(input,flux_input,fe_var.active_var_names());

      return this->build_neumman_func_common( input, system, fe_var, flux_input );
    }

    virtual std::string flux_input(const std::string& section ) const
    { return section+"/"+BoundaryConditionNames::bc_flux_var(); }

  };

  template<typename FunctionType>
  class ParsedTractionBCFactory : public ParsedFunctionNeumannBCFactory<FunctionType>
  {
  public:
    ParsedTractionBCFactory( const std::string& bc_type_name )
      : ParsedFunctionNeumannBCFactory<FunctionType>(bc_type_name)
    {}

    ~ParsedTractionBCFactory(){};

  protected:

    virtual std::string flux_input(const std::string& section ) const
    { return section+"/"+BoundaryConditionNames::traction_var(); }

  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_NEUMANN_BC_FACTORY_H
