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

#ifndef GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H
#define GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_old_style_base.h"
#include "grins/string_utils.h"

// libMesh
#include "libmesh/const_function.h"

namespace GRINS
{
  class IsothermalDirichletOldStyleBCFactory : public DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    IsothermalDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}

    ~IsothermalDirichletOldStyleBCFactory(){};

  protected:

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

  };

  std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
  inline
  IsothermalDirichletOldStyleBCFactory::build_func( const GetPot& input,
                                                    MultiphysicsSystem& /*system*/,
                                                    std::vector<std::string>& var_names,
                                                    const std::string& section )
  {
    libmesh_assert_equal_to(DirichletBCFactoryAbstract::_bc_ids->size(), 1 );
    libmesh_assert_equal_to(var_names.size(), 1 );

    std::string bc_id_string = StringUtilities::T_to_string<BoundaryID>( *(_bc_ids->begin()) );

    std::string input_var = section+"/T_wall_"+bc_id_string;

    if( !input.have_variable(input_var) )
      libmesh_error_msg("ERROR: Could not find input variable "+input_var+"!");

    libMesh::Number value = input(input_var, 0.0);

    return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( new libMesh::ConstFunction<libMesh::Number>(value) );
  }

} // end namespace GRINS

#endif // GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H
