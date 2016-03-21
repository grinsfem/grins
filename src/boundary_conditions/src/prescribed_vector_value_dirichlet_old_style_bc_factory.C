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

// This class
#include "grins/prescribed_vector_value_dirichlet_old_style_bc_factory.h"
#include "grins/string_utils.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/const_function.h"

namespace GRINS
{
  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Number> >
  PrescribedVectorValueDirichletOldStyleBCFactory::build_func( const GetPot& input,
                                                               MultiphysicsSystem& system,
                                                               std::vector<std::string>& var_names,
                                                               const std::string& section )
  {
    libmesh_assert_equal_to(DirichletBCFactoryAbstract::_bc_ids->size(), 1 );

    std::string bc_id_string = StringUtilities::T_to_string<BoundaryID>( *(_bc_ids->begin()) );

    std::string var_input_string = this->var_input_string();

    std::string input_string = section+"/"+var_input_string+"_"+bc_id_string;

    unsigned int n_comps = input.vector_variable_size(input_string);

    if( var_names.size() > n_comps )
      libmesh_error_msg("ERROR: Insufficient number of variable components in "+input_string+"!");

    libMesh::UniquePtr<libMesh::CompositeFunction<libMesh::Number> >
      remapped_func( new libMesh::CompositeFunction<libMesh::Number> );

    this->add_funcs(input,system,input_string,var_names,*(remapped_func.get()));

    return libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Number> >(remapped_func.release());
  }

  void PrescribedVectorValueDirichletOldStyleBCFactory::add_funcs( const GetPot& input,
                                                                   MultiphysicsSystem& system,
                                                                   const std::string& input_string,
                                                                   const std::vector<std::string>& var_names,
                                                                   libMesh::CompositeFunction<libMesh::Number>& composite_func ) const
  {
    for( unsigned int n = 0; n < var_names.size(); n++ )
      {
        std::vector<VariableIndex> dbc_vars(1,system.variable_number(var_names[n]));
        libMesh::Number value = input(input_string, 0.0, n);
        libMesh::ConstFunction<libMesh::Number> const_func(value);
        composite_func.attach_subfunction(const_func, dbc_vars);
      }
  }
  // Instantiate factories
  PrescribedVelDirichletOldStyleBCFactory grins_factory_prescribed_vel_old_style("prescribed_vel_old_style");
  PrescribedDispDirichletOldStyleBCFactory grins_factory_constant_displacement_old_style("constant_displacement_old_style");
  PrescribedSpeciesDirichletOldStyleBCFactory grins_factory_prescribed_species_old_style("prescribed_species_old_style");
} // end namespace GRINS
