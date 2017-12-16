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
#include "grins/velocity_penalty_base.h"

// GRINS
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{

  template<class Mu>
  VelocityPenaltyBase<Mu>::VelocityPenaltyBase( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    base_physics_name("VelocityPenalty"),
    _quadratic_scaling(false),
    _input(input)
  {
    if (this->_physics_name == PhysicsNaming::velocity_penalty2() ||
        this->_physics_name == PhysicsNaming::velocity_penalty2_adjoint_stab())
      base_physics_name.push_back('2');

    if (this->_physics_name == PhysicsNaming::velocity_penalty3() ||
        this->_physics_name == PhysicsNaming::velocity_penalty3_adjoint_stab())
      base_physics_name.push_back('3');

    this->read_input_options(input);
  }

  template<class Mu>
  void VelocityPenaltyBase<Mu>::read_input_options( const GetPot& input )
  {
    _quadratic_scaling =
      input("Physics/"+base_physics_name+"/quadratic_scaling", false);
  }

  template<class Mu>
  void VelocityPenaltyBase<Mu>::set_time_evolving_vars ( libMesh::FEMSystem* system)
  {
    libMesh::ParsedFEMFunction<libMesh::Number> *nvf
      (new libMesh::ParsedFEMFunction<libMesh::Number> (*system, ""));
    this->normal_vector_function.reset(nvf);

    this->set_parameter(*nvf, _input,
                        "Physics/"+base_physics_name+"/penalty_function",
                        this->zero_vector_function);

    libMesh::ParsedFEMFunction<libMesh::Number> *bvf
      (new libMesh::ParsedFEMFunction<libMesh::Number> (*system, ""));
    this->base_velocity_function.reset(bvf);

    this->set_parameter(*bvf, _input,
                        "Physics/"+base_physics_name+"/base_velocity",
                        this->zero_vector_function);

    if ((nvf->expression() == this->zero_vector_function) &&
        (bvf->expression() == this->zero_vector_function))
      std::cout << "Warning! Zero VelocityPenalty specified!" << std::endl;
  }


  template<class Mu>
  bool VelocityPenaltyBase<Mu>::compute_force
  ( const libMesh::Point& point,
    const AssemblyContext& context,
    const libMesh::NumberVectorValue& U,
    libMesh::NumberVectorValue& F,
    libMesh::NumberTensorValue *dFdU)
  {
    // Velocity discrepancy (current velocity minus base velocity)
    // normal to constraint plane, scaled by constraint penalty
    // value
    libMesh::DenseVector<libMesh::Number> output_vec(3);

    (*this->normal_vector_function)(context, point, context.time, output_vec);

    libMesh::NumberVectorValue U_N(output_vec(0),
                                   output_vec(1),
                                   output_vec(2));

    (*this->base_velocity_function)(context, point, context.time, output_vec);

    const libMesh::NumberVectorValue U_B(output_vec(0),
                                         output_vec(1),
                                         output_vec(2));

    const libMesh::NumberVectorValue U_Rel = U-U_B;

    // Old code
    // const libMesh::NumberVectorValue F1 = (U_Rel*U_N)*U_N; //

    // With correct sign and more natural normalization
    const libMesh::Number U_N_mag = std::sqrt(U_N*U_N);

    if (!U_N_mag)
      return false;

    const libMesh::NumberVectorValue U_N_unit = U_N/U_N_mag;

    F = -(U_Rel*U_N)*U_N_unit;

    if (dFdU)
      for (unsigned int i=0; i != 3; ++i)
        for (unsigned int j=0; j != 3; ++j)
          (*dFdU)(i,j) = -(U_N(j))*U_N_unit(i);

    // With quadratic scaling
    if (_quadratic_scaling)
      {
        const libMesh::Number U_Rel_mag = std::sqrt(U_Rel * U_Rel);

        // Modify dFdU first so as to reuse the old value of F
        if (dFdU)
          {
            // dU_Rel/dU = I
            // d(U_Rel*U_Rel)/dU = 2*U_Rel
            // d|U_Rel|/dU = U_Rel/|U_Rel|

            const libMesh::NumberVectorValue U_Rel_unit = U_Rel/U_Rel_mag;

            (*dFdU) *= U_Rel_mag;

            if (U_Rel_mag)
              for (unsigned int i=0; i != 3; ++i)
                for (unsigned int j=0; j != 3; ++j)
                  (*dFdU)(i,j) += F(i)*U_Rel_unit(j);
          }

        F *= U_Rel_mag;
      }

    // With correction term to avoid doing work on flow
    bool do_correction = false;
    if (do_correction)
      {
        if (dFdU)
          libmesh_not_implemented();

        const libMesh::Number U_Rel_mag_sq = U_Rel * U_Rel;
        if (U_Rel_mag_sq)
          {
            F -= (U_Rel*F)*U_Rel/U_Rel_mag_sq;
          }
      }
    return true;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(VelocityPenaltyBase);
