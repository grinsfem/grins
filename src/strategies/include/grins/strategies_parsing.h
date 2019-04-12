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

#ifndef GRINS_STRATEGIES_PARSING_H
#define GRINS_STRATEGIES_PARSING_H

// C++
#include <string>

// Forward declarations
class GetPot;
namespace libMesh
{
  class SystemNorm;
}
namespace GRINS
{
  class StrategiesParsing
  {
  public:

    StrategiesParsing(){};

    ~StrategiesParsing(){};

    static int extra_quadrature_order( const GetPot& input );

    static std::string adjoint_residual_error_estimator()
    { return "adjoint_residual"; }

    static std::string adjoint_refinement_error_estimator()
    { return "adjoint_refinement"; }

    static std::string kelly_error_estimator()
    { return "kelly"; }

    static std::string patch_recovery_error_estimator()
    { return "patch_recovery"; }

    //! Option to let user manually trigger adjoint solve
    static bool do_adjoint_solve( const GetPot& input );
  };

} // end namespace GRINS

#endif // GRINS_STRATEGIES_PARSING_H
