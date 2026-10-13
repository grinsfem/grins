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

#ifndef GRINS_PENALTY_DIRICHLET_H
#define GRINS_PENALTY_DIRICHLET_H

// GRINS
#include "grins/boundary_restricted.h"
#include "grins/physics.h"

namespace GRINS
{

  class PenaltyDirichlet : public Physics, public BoundaryRestricted
  {
  public:

    PenaltyDirichlet( const PhysicsName& physics_name, const GetPot& input );

    virtual ~PenaltyDirichlet() = default;

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext & context ) override;

    //! Initialize pinning helper object
    virtual void auxiliary_init( MultiphysicsSystem & system ) override;

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time independent part(s) of physics for element sides
    virtual void side_constraint( bool compute_jacobian,
                                  AssemblyContext & context ) override;

  protected:

    std::string _variablename_to_pin;

    unsigned int _variable_to_pin;

    libMesh::Real _penalty;

    libMesh::Real _target;
  };

} // namespace GRINS

#endif // GRINS_PENALTY_DIRICHLET_H
