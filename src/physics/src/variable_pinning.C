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
#include "grins/variable_pinning.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"

namespace GRINS
{

  VariablePinning::VariablePinning( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _var_pinning(input,physics_name),
      _penalty(input("Physics/" + PhysicsNaming::variable_pinning() + "/penalty", 1e10)),
      _pin_variable(false)
  {
    std::string pin_variable_str = "Physics/" +
      PhysicsNaming::variable_pinning() + "/pin_variable";
    if (input.have_variable(pin_variable_str))
      {
        _variablename_to_pin = input(pin_variable_str, std::string());
        _pin_variable = true;
      }
  }

  void VariablePinning::auxiliary_init( MultiphysicsSystem & system )
  {
    if( _pin_variable )
      {
        _variable_to_pin =
          system.variable_number(_variablename_to_pin);
        _var_pinning.check_pin_location(system.get_mesh());
      }
    else
      {
        // Commenting out pin_variable makes testing easier but never
        // makes sense in a real run.
        libMesh::out <<
          "Warning!  VariablePinning physics requested but not used!"
                     << std::endl;
      }
  }

  void VariablePinning::init_context( AssemblyContext & context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    if ( _pin_variable )
      {
        context.get_element_fe(_variable_to_pin)->get_JxW();
        context.get_element_fe(_variable_to_pin)->get_phi();
      }
  }

  void VariablePinning::element_constraint
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // Pin var = var_value at p_point
    if( _pin_variable )
      _var_pinning.pin_value( context, compute_jacobian, this->_variable_to_pin, this->_penalty );
  }

} // namespace GRINS
