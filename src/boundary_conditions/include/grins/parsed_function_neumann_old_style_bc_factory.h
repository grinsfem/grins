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

#ifndef GRINS_PARSED_FUNCTION_NEUMANN_OLD_STYLE_BC_FACTORY_H
#define GRINS_PARSED_FUNCTION_NEUMANN_OLD_STYLE_BC_FACTORY_H

// GRINS
#include "grins/neumann_bc_old_style_factory_abstract.h"
#include "grins/parsed_function_neumann_bc_factory_helper.h"
#include "grins/boundary_condition_names.h"
#include "grins/string_utils.h"

namespace GRINS
{
  template<typename FunctionType>
  class ParsedFunctionNeumannOldStyleBCFactory : public NeumannBCOldStyleFactoryAbstract,
                                                 public ParsedFunctionNeumannBCFactoryHelper<FunctionType>
  {
  public:

    ParsedFunctionNeumannOldStyleBCFactory( const std::string& bc_type_name )
      : NeumannBCOldStyleFactoryAbstract(bc_type_name),
        ParsedFunctionNeumannBCFactoryHelper<FunctionType>()
    {}

    ~ParsedFunctionNeumannOldStyleBCFactory(){};

  protected:

    virtual std::shared_ptr<NeumannBCAbstract>
    build_neumann_func( const GetPot& input,
                        MultiphysicsSystem& system,
                        const FEVariablesBase& fe_var,
                        const std::string& section );

    virtual std::string flux_input() const =0;

  };

  template<typename FunctionType>
  inline
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannOldStyleBCFactory<FunctionType>::build_neumann_func( const GetPot& input,
                                                                            MultiphysicsSystem& system,
                                                                            const FEVariablesBase& fe_var,
                                                                            const std::string& section )
  {
    libmesh_assert_equal_to( this->_bc_ids->size(), 1 );

    std::string flux_input = section+"/"+this->flux_input()+"_"+
      StringUtilities::T_to_string<unsigned int>( *(this->_bc_ids->begin()) );

    // Make sure flux input specified and consistent with var_names size
    this->check_for_flux(input,flux_input,fe_var.active_var_names());

    return this->build_neumman_func_common( input, system, fe_var, flux_input );
  }

  template<typename FunctionType>
  class TractionOldStyleBCFactory : public ParsedFunctionNeumannOldStyleBCFactory<FunctionType>
  {
  public:

    TractionOldStyleBCFactory( const std::string& bc_type_name )
      : ParsedFunctionNeumannOldStyleBCFactory<FunctionType>(bc_type_name)
    {}

    ~TractionOldStyleBCFactory(){};

  protected:

    virtual std::string flux_input() const
    { return "traction"; }
  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_NEUMANN_OLD_STYLE_BC_FACTORY_H
