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

// These functions
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/system_norm.h"

namespace GRINS
{

  int StrategiesParsing::extra_quadrature_order( const GetPot& input )
  {
    int extra_order = input("Strategies/Assembly/extra_quadrature_order", 0);

    if( extra_order < 0 )
      libmesh_error_msg("ERROR: extra_quadrature_order must be non-negative!");

    return extra_order;
  }

  bool StrategiesParsing::do_adjoint_solve( const GetPot& input )
  {
    bool do_adjoint_solve = false;
    std::string old_option = "linear-nonlinear-solver/do_adjoint_solve";
    std::string new_option = "Strategies/Adjoint/do_adjoint_solve";

    bool have_old_option = input.have_variable(old_option);
    bool have_new_option = input.have_variable(new_option);
    if( have_old_option  &&  have_new_option )
      libmesh_error_msg("ERROR: Cannot specify both "+old_option+" and "+new_option+"!");

    if( have_old_option )
      do_adjoint_solve = input( old_option, false );

    if( have_new_option )
      do_adjoint_solve = input( new_option, false );

    return do_adjoint_solve;
  }
} // end namespace GRINS
