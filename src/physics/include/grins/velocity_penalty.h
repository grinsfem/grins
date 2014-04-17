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


#ifndef GRINS_VELOCITY_PENALTY_H
#define GRINS_VELOCITY_PENALTY_H

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/cached_values.h"
#include "grins/inc_navier_stokes_base.h"

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/getpot.h"

// C++
#include <string>

namespace GRINS
{

  //! Physics class for Velocity Penalty
  /*
    This physics class imposes a penalty on any velocity component in
    the direction of (and proportional to) a specified vector field.
   */
  class VelocityPenalty : public IncompressibleNavierStokesBase
  {
  public:

    VelocityPenalty( const std::string& physics_name, const GetPot& input );

    ~VelocityPenalty();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );
    
    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Constraint part(s)
    virtual void element_constraint( bool compute_jacobian,
				     AssemblyContext& context,
				     CachedValues& cache );

  private:

    AutoPtr<libMesh::FunctionBase<Number> > normal_vector_function;

    AutoPtr<libMesh::FunctionBase<Number> > base_velocity_function;

    VelocityPenalty();
  };

} // end namespace block

#endif // GRINS_VELOCITY_PENALTY_H
