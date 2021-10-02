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


// This class
#include "grins/parsed_qoi_base.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/parsed_fem_function.h"

// C++
#include <string>

namespace GRINS
{
  ParsedQoIBase::ParsedQoIBase(const ParsedQoIBase & original)
    : QoIBase(original)
  {
    if (original.qoi_functional.get())
      {
        this->qoi_functional = original.qoi_functional->clone();

        this->move_parameter
          ( *libMesh::cast_ptr<libMesh::ParsedFEMFunction<libMesh::Number>*>
            (original.qoi_functional.get()),
            *libMesh::cast_ptr<libMesh::ParsedFEMFunction<libMesh::Number>*>
            (this->qoi_functional.get()) );
      }
  }

  void ParsedQoIBase::init_qoi_functional( const GetPot & input,
                                           const MultiphysicsSystem & system,
                                           const std::string & input_string )
  {
    std::unique_ptr<libMesh::ParsedFEMFunction<libMesh::Number>> qf =
      libmesh_make_unique<libMesh::ParsedFEMFunction<libMesh::Number>>
       (system, "");

    this->set_parameter(*qf, input,
                        input_string, "DIE!");

    this->get_var_indices( qf->expression(), system, _var_indices );

    this->qoi_functional = std::move(qf);
  }

  void ParsedQoIBase::get_var_indices( const std::string & expression,
                                       const MultiphysicsSystem & system,
                                       std::set<unsigned int> & var_indices )
  {
    std::vector<unsigned int> all_var_indices;
    system.get_all_variable_numbers(all_var_indices);

    for( auto var : all_var_indices )
      {
        const std::string & var_name = system.variable_name(var);
        if( expression.find(var_name) != std::string::npos )
          var_indices.insert(var);
      }
  }

} //namespace GRINS
