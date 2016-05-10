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

// These classes
#include "grins/parsed_function_neumann_old_style_bc_factory.h"

namespace GRINS
{
  template<typename FunctionType>
  SharedPtr<NeumannBCAbstract>
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

  template class ParsedFunctionNeumannOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >;
  template class ParsedFunctionNeumannOldStyleBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >;

  TractionOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >
  grins_factory_traction_old_style_functionbase("constant_traction_old_style");

} // end namespace GRINS
