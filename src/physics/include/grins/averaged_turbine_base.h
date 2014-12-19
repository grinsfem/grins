//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "libmesh/function_base.h"
#include "libmesh/getpot.h"
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

    ~AveragedTurbineBase();

    //! Initialization of variables
    /*!
      Add turbine_speed variable to system; call base function to
      initialize Navier-Stokes variables.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets turbine_speed and velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    bool compute_force ( const libMesh::Point& point,
                         const libMesh::Real time,
                         const libMesh::NumberVectorValue& U,
                         libMesh::Number s,
                         libMesh::NumberVectorValue& U_B_1,
                         libMesh::NumberVectorValue& F,
                         libMesh::NumberTensorValue *dFdU = NULL,
                         libMesh::NumberVectorValue *dFds = NULL);
 
    VariableIndex fan_speed_var() const { return _fan_speed_var; }

  protected:

    // ``Base'' velocity of the moving fan blades as a function of x,y,z
    //
    // This velocity will be scaled by the fan_speed variable
    // during the simulation; values in base_velocity_function should
    // therefore correspond to a fan speed of 1 radian per second.
    // Iff there are no fan blades moving past a location, the
    // base velocity there is specified as zero.
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > base_velocity_function;

    // "Up" direction of fan airflow, as a function of x,y,z
    // For most fans this will be a constant, the axis of rotation.
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > local_vertical_function;

    // Coefficients of lift and drag as a function of angle "t" in
    // radians.  Should be well defined on [-pi, pi].
    //
    // No, "t" is not time in these functions.
    // Yes, I'm abusing FunctionBase.
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > lift_function;
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > drag_function;

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
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > torque_function;

    // Moment of inertia of the spinning component of the turbine
    // (measured in kg-m^2)
    libMesh::Number moment_of_inertia;

    // Initial speed of the spinning component of the turbine
    // (measured in rad/s)
    libMesh::Number initial_speed;

    // The chord length of the fan wing cross-section.  For fan blades
    // with constant cross-section this will be a constant.
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > chord_function;

    // The area swept out by the local fan wing cross-section.  For
    // cylindrical areas swept out by N fan blades, this is just
    // pi*r^2*h/N
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > area_swept_function;

    // The angle-of-attack between the fan wing chord line and the fan
    // velocity vector, in radians.  For fan blades with no "twist"
    // this will be a constant.
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > aoa_function;

    VariableIndex _fan_speed_var; /* Index for turbine speed scalar */

    std::string _fan_speed_var_name;

  private:

    AveragedTurbineBase();
  };

} // end namespace block

#endif // GRINS_AVERAGED_TURBINE_BASE_H
