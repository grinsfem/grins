//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/velocity_drag_base.h"

// GRINS
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/parsed_function.h"

namespace GRINS
{

  template<class Mu>
  VelocityDragBase<Mu>::VelocityDragBase( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    _coefficient("")
  {
    this->read_input_options(input);
  }

  template<class Mu>
  void VelocityDragBase<Mu>::read_input_options( const GetPot& input )
  {
    this->set_parameter
      (_exponent, input, "Physics/"+PhysicsNaming::velocity_drag()+"/exponent", 2);

    this->set_parameter(_coefficient, input,
                        "Physics/"+PhysicsNaming::velocity_drag()+"/coefficient",
                        "0");

    if (_coefficient.expression() == "0")
      libmesh_error_msg("Warning! Zero VelocityDrag specified!" <<
                        std::endl);
  }

  template<class Mu>
  bool VelocityDragBase<Mu>::compute_force
  ( const libMesh::Point& point,
    const libMesh::Real time,
    const libMesh::NumberVectorValue& U,
    libMesh::NumberVectorValue& F,
    libMesh::NumberTensorValue *dFdU)
  {
    const libMesh::Number Umag = U.norm();

    const libMesh::Number coeff_val = _coefficient(point, time);

    if (coeff_val == 0)
      return false;

    const libMesh::Number F_coeff = std::pow(Umag, _exponent-1) * -coeff_val;

    F = F_coeff * U;

    if (dFdU)
      {
        const libMesh::Number J_coeff =
          std::pow(Umag, _exponent-2) * -coeff_val * (_exponent-1);

        const libMesh::Number invUmag = 1/Umag;

        const libMesh::Number JinvU = J_coeff*invUmag;

        for (unsigned int i=0; i != 3; ++i)
          {
            const libMesh::Number JinvUI = JinvU*U(i);
            for (unsigned int j=0; j != 3; ++j)
              (*dFdU)(i,j) = JinvUI*U(j);

            (*dFdU)(i,i) += F_coeff;
          }
      }

    return true;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(VelocityDragBase);
