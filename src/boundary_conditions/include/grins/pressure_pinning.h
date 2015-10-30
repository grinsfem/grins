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
  class MeshBase;
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

    //! Check the mesh to ensure pin location is found
    /*! If the pin location is found, we set _pin_location_found to true
        (for assertion later). If not found, we throw an error. */
    void check_pin_location( const libMesh::MeshBase& mesh );

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

    //! Cache element id for element that contains _pin_location
    /*! We will initalize this to libMesh::DofObject::invalid_id
        and use that to check whether or not we located an element
        that contains the _pin_location. */
    libMesh::dof_id_type _pinned_elem_id;

  };
}
#endif //PRESSURE_PINNING_H
