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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/ic_handling_base.h"

// GRINS
#include "grins/string_utils.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/dof_map.h"
#include "libmesh/parsed_function.h"
#include "libmesh/const_function.h"

// C++
#include "sstream"

namespace GRINS
{
  ICHandlingBase::ICHandlingBase(const std::string& physics_name)
    : _physics_name( physics_name )
  {
    return;
  }

  ICHandlingBase::~ICHandlingBase()
  {
    return;
  }

  void ICHandlingBase::attach_initial_func
  ( const libMesh::FunctionBase<libMesh::Number>& initial_val)
  {
    _ic_func = initial_val.clone();
  }

  void ICHandlingBase::read_ic_data( const GetPot& input, const std::string& id_str,
                                     const std::string& ic_str,
                                     const std::string& var_str,
                                     const std::string& value_str)
  {
    int num_ids = input.vector_variable_size(id_str);
    int num_ics = input.vector_variable_size(ic_str);
    int num_vars = input.vector_variable_size(var_str);
    int num_values = input.vector_variable_size(value_str);

    if( num_ids != num_ics )
      {
        std::cerr << "Error: number of subdomain ids " << num_ids
                  << "must equal number of initial condition types " << num_ics
                  << std::endl;
        libmesh_error();
      }

    if( num_ids != num_vars )
      {
        std::cerr << "Error: number of subdomain ids " << num_ids
                  << "must equal number of variable name lists " << num_vars
                  << std::endl;
        libmesh_error();
      }

    if( num_ids != num_values )
      {
        std::cerr << "Error: number of subdomain ids " << num_ids
                  << "must equal number of initial condition values " << num_values
                  << std::endl;
        libmesh_error();
      }

    if( num_ids > 1 )
      {
        std::cerr << "Error: GRINS does not yet support per-subdomain initial conditions" << std::endl;
        libmesh_not_implemented();
      }

    for( int i = 0; i < num_ids; i++ )
      {
        int ic_id = input(id_str, -1, i );
        std::string ic_type_in = input(ic_str, "NULL", i );
        std::string ic_value_in = input(value_str, "NULL", i );
        std::string ic_vars_in = input(var_str, "NULL", i );

        int ic_type = this->string_to_int( ic_type_in );

        std::stringstream ss;
        ss << ic_id;
        std::string ic_id_string = ss.str();

        this->init_ic_types( ic_id, ic_id_string, ic_type, ic_vars_in, ic_value_in, input );
      }

    return;
  }

  void ICHandlingBase::init_ic_data( const libMesh::FEMSystem& system,
                                     libMesh::CompositeFunction<libMesh::Number>& all_ics )
  {
    if (this->get_ic_func())
      {
        std::vector<unsigned int> index_map;

        libmesh_assert(_subfunction_variables.size());

        for (unsigned int i=0; i != _subfunction_variables.size();
             ++i)
          index_map.push_back
            (system.variable_number(_subfunction_variables[i]));

        all_ics.attach_subfunction(*this->get_ic_func(), index_map);
      }
  }

  int ICHandlingBase::string_to_int( const std::string& ic_type_in ) const
  {
    int ic_type_out;
    if( ic_type_in == "parsed" )
      {
        ic_type_out = PARSED;
      }
    else if( ic_type_in == "constant" )
      {
        ic_type_out = CONSTANT;
      }
    else
      {
        std::cerr << "=========================================================="  << std::endl
                  << "Error: Invalid ic_type " << ic_type_in                       << std::endl
                  << "       Physics class is " << _physics_name                   << std::endl
                  << "=========================================================="  << std::endl;
        libmesh_error();
      }
    return ic_type_out;
  }

  void ICHandlingBase::init_ic_types( const libMesh::subdomain_id_type /*ic_id*/,
                                      const std::string& /*ic_id_string*/,
                                      const int ic_type,
                                      const std::string& ic_vars_string,
                                      const std::string& ic_value_string,
                                      const GetPot& /*input*/ )
  {
    StringUtilities::split_string(ic_vars_string, ":", _subfunction_variables);

    libmesh_assert(_subfunction_variables.size());

    switch(ic_type)
      {
      case(PARSED):
        {
          _ic_func = std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
            (new libMesh::ParsedFunction<libMesh::Number>(ic_value_string));
        }
        break;

      case(CONSTANT):
        {
          _ic_func = std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
            (new libMesh::ConstFunction<libMesh::Number>
             (StringUtilities::string_to_T<libMesh::Number>(ic_value_string)));
        }
        break;

      default:
        {
          std::cerr << "=========================================================="
                    << "Error: Invalid IC type for " << _physics_name << std::endl
                    << "       Detected IC type was " << ic_type << std::endl
                    << "==========================================================" << std::endl;
          libmesh_error();
        }
      }
    return;
  }

} // namespace GRINS
