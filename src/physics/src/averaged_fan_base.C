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
#include "grins/averaged_fan_base.h"

// GRINS
#include "grins/inc_nav_stokes_macro.h"

// libMesh

namespace GRINS
{

  template<class Mu>
  AveragedFanBase<Mu>::AveragedFanBase( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    base_velocity_function(""),
    local_vertical_function(""),
    lift_function(""),
    drag_function(""),
    chord_function(""),
    area_swept_function(""),
    aoa_function("")
  {
    this->read_input_options(input);
  }

  template<class Mu>
  void AveragedFanBase<Mu>::read_input_options( const GetPot& input )
  {
    this->set_parameter(base_velocity_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/base_velocity",
                        this->zero_vector_function);

    if (base_velocity_function.expression() == this->zero_vector_function)
      libmesh_error_msg("Error! Zero AveragedFan specified!" <<
                        std::endl);

    this->set_parameter(local_vertical_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/local_vertical",
                        this->zero_vector_function);

    if (local_vertical_function.expression() == this->zero_vector_function)
      libmesh_error_msg("Error! Zero LocalVertical specified!" <<
                        std::endl);

    this->set_parameter(lift_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/lift",
                        "0");

    if (lift_function.expression() == "0")
      std::cout << "Warning! Zero lift function specified!" << std::endl;

    this->set_parameter(drag_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/drag",
                        "0");

    if (drag_function.expression() == "0")
      std::cout << "Warning! Zero drag function specified!" << std::endl;

    this->set_parameter(chord_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/chord_length",
                        "0");

    if (chord_function.expression() == "0")
      libmesh_error_msg("Error! Zero chord function specified!" <<
                        std::endl);

    this->set_parameter(area_swept_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/area_swept",
                        "0");

    if (area_swept_function.expression() == "0")
      libmesh_error_msg("Error! Zero area_swept_function specified!" <<
                        std::endl);

    this->set_parameter(aoa_function, input,
                        "Physics/"+PhysicsNaming::averaged_fan()+"/angle_of_attack",
                        "00000");

    if (aoa_function.expression() == "00000")
      libmesh_error_msg("Error! No angle-of-attack specified!" <<
                        std::endl);
  }

  template<class Mu>
  bool AveragedFanBase<Mu>::compute_force
  ( const libMesh::Point& point,
    const libMesh::Real time,
    const libMesh::NumberVectorValue& U,
    libMesh::NumberVectorValue& F,
    libMesh::NumberTensorValue *dFdU)
  {
    // Find base velocity of moving fan at this point
    libMesh::DenseVector<libMesh::Number> output_vec(3);

    base_velocity_function(point, time,
                           output_vec);

    const libMesh::NumberVectorValue U_B(output_vec(0),
                                         output_vec(1),
                                         output_vec(2));

    const libMesh::Number U_B_size = U_B.norm();

    // If there's no base velocity there's no fan
    if (!U_B_size)
      return false;

    // Normal in fan velocity direction
    const libMesh::NumberVectorValue N_B =
      libMesh::NumberVectorValue(U_B/U_B_size);

    local_vertical_function(point, time,
                            output_vec);

    // Normal in fan vertical direction
    const libMesh::NumberVectorValue N_V(output_vec(0),
                                         output_vec(1),
                                         output_vec(2));

    // Normal in radial direction (or opposite radial direction,
    // for fans turning clockwise!)
    const libMesh::NumberVectorValue N_R = N_B.cross(N_V);

    // Fan-wing-plane component of local relative velocity
    const libMesh::NumberVectorValue U_P = U - (U*N_R)*N_R - U_B;

    const libMesh::Number U_P_size = U_P.norm();

    // If there's no flow in the fan's frame of reference, there's no
    // lift or drag.  FIXME - should we account for drag in the
    // out-of-plane direction?
    if (!U_P_size)
      return false;

    // Direction opposing drag
    const libMesh::NumberVectorValue N_drag =
      libMesh::NumberVectorValue(-U_P/U_P_size);

    // Direction opposing lift
    const libMesh::NumberVectorValue N_lift = N_drag.cross(N_R);

    // "Forward" velocity
    const libMesh::Number u_fwd = -(U_P * N_B);

    // "Upward" velocity
    const libMesh::Number u_up = U_P * N_V;

    // If there's no forward or upward velocity we should have already
    // returned false
    libmesh_assert (u_up || u_fwd);

    // Angle WRT fan velocity direction
    const libMesh::Number part_angle = std::atan2(u_up, u_fwd);

    // Angle WRT fan chord
    const libMesh::Number angle = part_angle +
      aoa_function(point, time);

    const libMesh::Number C_lift  = lift_function(point, angle);
    const libMesh::Number C_drag  = drag_function(point, angle);

    const libMesh::Number chord = chord_function(point, time);
    const libMesh::Number area  = area_swept_function(point, time);

    const libMesh::Number v_sq = U_P*U_P;

    const libMesh::Number LDfactor = 0.5 * this->_rho * v_sq * chord / area;
    const libMesh::Number lift = C_lift * LDfactor;
    const libMesh::Number drag = C_drag * LDfactor;

    // Force
    F = lift * N_lift + drag * N_drag;

    if (dFdU)
      {
        // FIXME: Jacobians here are very inexact!
        // Dropping all AoA dependence on U terms!
        const libMesh::NumberVectorValue LDderivfactor =
          (N_lift*C_lift+N_drag*C_drag) *
          this->_rho * chord / area;

        for (unsigned int i=0; i != 3; ++i)
          for (unsigned int j=0; j != 3; ++j)
            (*dFdU)(i,j) = LDderivfactor(i) * U_P(j);
      }

    return true;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(AveragedFanBase);
