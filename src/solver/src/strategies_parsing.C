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

// These functions
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/system_norm.h"

namespace GRINS
{
  bool StrategiesParsing::is_mesh_adaptive( const GetPot& input )
  {
    return input("MeshAdaptivity/mesh_adaptive", false );
  }

  double StrategiesParsing::parse_target_tolerance( const GetPot& input )
  {
    return input("unsteady-solver/target_tolerance", 0.0 );
  }

  double StrategiesParsing::parse_upper_tolerance( const GetPot& input )
  {
    return input("unsteady-solver/upper_tolerance", 0.0 );
  }

  double StrategiesParsing::parse_max_growth( const GetPot& input )
  {
    return input("unsteady-solver/max_growth", 0.0 );
  }

  void StrategiesParsing::parse_component_norm( const GetPot& input, libMesh::SystemNorm& component_norm )
  {
    const unsigned int n_component_norm =
      input.vector_variable_size("unsteady-solver/component_norm");
    for (unsigned int i=0; i != n_component_norm; ++i)
      {
        const std::string current_norm = input("component_norm", std::string("L2"), i);
        // TODO: replace this with string_to_enum with newer libMesh
        if (current_norm == "GRINSEnums::L2")
          component_norm.set_type(i, libMesh::L2);
        else if (current_norm == "GRINSEnums::H1")
          component_norm.set_type(i, libMesh::H1);
        else
          libmesh_not_implemented();
      }
  }

  int StrategiesParsing::extra_quadrature_order( const GetPot& input )
  {
    int extra_order = input("Strategies/Assembly/extra_quadrature_order", 0);

    if( extra_order < 0 )
      libmesh_error_msg("ERROR: extra_quadrature_order must be non-negative!");

    return extra_order;
  }

} // end namespace GRINS
