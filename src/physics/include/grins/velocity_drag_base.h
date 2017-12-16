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


#ifndef GRINS_VELOCITY_DRAG_BASE_H
#define GRINS_VELOCITY_DRAG_BASE_H

// GRINS
#include "grins_config.h"
#include "grins/inc_navier_stokes_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parsed_function.h"
#include "libmesh/tensor_value.h"

// C++
#include <string>

namespace GRINS
{

  template<class Viscosity>
  class VelocityDragBase : public IncompressibleNavierStokesBase<Viscosity>
  {
  public:

    VelocityDragBase( const std::string& physics_name, const GetPot& input );

    ~VelocityDragBase(){};

    bool compute_force ( const libMesh::Point& point,
                         const libMesh::Real time,
                         const libMesh::NumberVectorValue& U,
                         libMesh::NumberVectorValue& F,
                         libMesh::NumberTensorValue *dFdU = NULL);

  protected:

    libMesh::Real _exponent;
    libMesh::ParsedFunction<libMesh::Number> _coefficient;

  private:

    VelocityDragBase();

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );
  };

} // end namespace block

#endif // GRINS_VELOCITY_DRAG_BASE_H
