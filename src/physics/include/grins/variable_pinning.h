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

#ifndef GRINS_VARIABLE_PINNING_H
#define GRINS_VARIABLE_PINNING_H

// GRINS
#include "grins/physics.h"
#include "grins/pressure_pinning.h"

namespace GRINS
{

  class VariablePinning : public Physics
  {

  public:

    VariablePinning( const GRINS::PhysicsName& physics_name, const GetPot& input );
    ~VariablePinning(){};

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext & context );

    //! Initialize pinning helper object
    virtual void auxiliary_init( MultiphysicsSystem & system );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

  protected:

    PressurePinning _var_pinning;

    std::string _variablename_to_pin;

    unsigned int _variable_to_pin;

    libMesh::Real _penalty;

    bool _pin_variable;

  private:

    VariablePinning();

  };

} // namespace GRINS

#endif // GRINS_VARIABLE_PINNING_H
