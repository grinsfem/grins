//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
#ifndef PRESSURE_PINNING_H
#define PRESSURE_PINNING_H

// libMesh stuff
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/point.h"

//GRINS
#include "grins/var_typedefs.h"

// libMesh forward declarations
class Getpot;

namespace libMesh
{
  class DiffContext;
}

namespace GRINS
{
  //! Class to hold typical boundary condition methods
  /*!
    This class holds functions to apply generic versions of
    Dirichlet and Neumann boundary conditions.
  */
  class PressurePinning
  {
  public:

    PressurePinning( const GetPot& input,
		     const std::string& physics_name );
    ~PressurePinning();

    /*! The idea here is to pin a variable to a particular value if there is
      a null space - e.g. pressure for IncompressibleNavierStokes. */
    void pin_value( libMesh::DiffContext& context,
		    const bool request_jacobian,
		    const GRINS::VariableIndex var,
		    const double penalty = 1.0 );

  private:

    //! Value of pressure we wish to pin
    libMesh::Number _pin_value;

    //! Location we want to pin the pressure
    libMesh::Point _pin_location;

  };
}
#endif //PRESSURE_PINNING_H
