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


#ifndef GRINS_AVERAGED_TURBINE_BASE_H
#define GRINS_AVERAGED_TURBINE_BASE_H

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
  class AveragedTurbineBase : public IncompressibleNavierStokesBase<Viscosity>
  {
  public:

    AveragedTurbineBase( const std::string& physics_name, const GetPot& input );

    virtual ~AveragedTurbineBase() = default;

    //! Sets turbine_speed and velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system ) override;

    bool compute_force ( const libMesh::Point& point,
                         const libMesh::Real time,
                         const libMesh::NumberVectorValue& U,
                         libMesh::Number s,
                         libMesh::NumberVectorValue& U_B_1,
                         libMesh::NumberVectorValue& F,
                         libMesh::NumberTensorValue *dFdU = NULL,
                         libMesh::NumberVectorValue *dFds = NULL);

    VariableIndex fan_speed_var() const { return _var.var(); }

  protected:

    // ``Base'' velocity of the moving fan blades as a function of x,y,z
    //
    // This velocity will be scaled by the fan_speed variable
    // during the simulation; values in base_velocity_function should
    // therefore correspond to a fan speed of 1 radian per second.
    // Iff there are no fan blades moving past a location, the
    // base velocity there is specified as zero.
    libMesh::ParsedFunction<libMesh::Number> base_velocity_function;

    // "Up" direction of fan airflow, as a function of x,y,z
    // For most fans this will be a constant, the axis of rotation.
    libMesh::ParsedFunction<libMesh::Number> local_vertical_function;

    // Coefficients of lift and drag as a function of angle "t" in
    // radians.  Should be well defined on [-pi, pi].
    //
    // No, "t" is not time in these functions.
    // Yes, I'm abusing FunctionBase.
    libMesh::ParsedFunction<libMesh::Number> lift_function;
    libMesh::ParsedFunction<libMesh::Number> drag_function;

    // Mechanical driving torque function (*including* non-fluid
    // friction losses!) on the turbine (signed, measured in
    // Newton-meters) as a function of angular velocity turbine speed
    // "t" (measured in rad/s).
    //
    // Should probably be carefully defined for t>0 and t<0 to avoid
    // potential unstable startup.
    //
    // Should theoretically be useable as power input (same sign as t)
    // to model propeller fans or output (opposite sign from t) to
    // model turbine fans.
    //
    // No, "t" is not time or angle of attack in this function.
    // Yes, I'm really abusing FunctionBase.
    libMesh::ParsedFunction<libMesh::Number> torque_function;

    // Moment of inertia of the spinning component of the turbine
    // (measured in kg-m^2)
    libMesh::Number moment_of_inertia;

    // Initial speed of the spinning component of the turbine
    // (measured in rad/s)
    libMesh::Number initial_speed;

    // The chord length of the fan wing cross-section.  For fan blades
    // with constant cross-section this will be a constant.
    libMesh::ParsedFunction<libMesh::Number> chord_function;

    // The area swept out by the local fan wing cross-section.  For
    // cylindrical areas swept out by N fan blades, this is just
    // pi*r^2*h/N
    libMesh::ParsedFunction<libMesh::Number> area_swept_function;

    // The angle-of-attack between the fan wing chord line and the fan
    // velocity vector, in radians.  For fan blades with no "twist"
    // this will be a constant.
    libMesh::ParsedFunction<libMesh::Number> aoa_function;

    ScalarVariable & _var;

  private:

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );

  };

} // end namespace block

#endif // GRINS_AVERAGED_TURBINE_BASE_H
