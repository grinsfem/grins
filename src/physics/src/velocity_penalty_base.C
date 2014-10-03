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


// This class
#include "grins/velocity_penalty_base.h"

// libMesh
#include "libmesh/parsed_function.h"
#include "libmesh/zero_function.h"

namespace GRINS
{

  VelocityPenaltyBase::VelocityPenaltyBase( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  VelocityPenaltyBase::~VelocityPenaltyBase()
  {
    return;
  }

  void VelocityPenaltyBase::read_input_options( const GetPot& input )
  {
    std::string penalty_function =
      input("Physics/"+velocity_penalty+"/penalty_function",
        std::string("0"));

    if (penalty_function == "0")
      this->normal_vector_function.reset
        (new libMesh::ZeroFunction<libMesh::Number>());
    else
      this->normal_vector_function.reset
        (new libMesh::ParsedFunction<libMesh::Number>(penalty_function));

    std::string base_function =
      input("Physics/"+velocity_penalty+"/base_velocity",
        std::string("0"));

    if (penalty_function == "0" && base_function == "0")
      std::cout << "Warning! Zero VelocityPenalty specified!" << std::endl;

    if (base_function == "0")
      this->base_velocity_function.reset
        (new libMesh::ZeroFunction<libMesh::Number>());
    else
      this->base_velocity_function.reset
        (new libMesh::ParsedFunction<libMesh::Number>(base_function));

  }

  bool VelocityPenaltyBase::compute_force
    ( const libMesh::Point& point,
      const libMesh::Real time,
      const libMesh::NumberVectorValue& U,
      libMesh::NumberVectorValue& F,
      libMesh::NumberTensorValue *dFdU)
  {
    // Velocity discrepancy (current velocity minus base velocity)
    // normal to constraint plane, scaled by constraint penalty
    // value
    libmesh_assert(normal_vector_function.get());
    libmesh_assert(base_velocity_function.get());

    libMesh::DenseVector<libMesh::Number> output_vec(3);

    (*normal_vector_function)(point, time,
                              output_vec);

    libMesh::NumberVectorValue U_N(output_vec(0),
                                   output_vec(1),
                                   output_vec(2));

    (*base_velocity_function)(point, time,
                              output_vec);

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

    // With correction term to avoid doing work on flow
    bool do_correction = false;
    if (do_correction)
      {
        const libMesh::Number U_Rel_mag_sq = U_Rel * U_Rel;
        if (U_Rel_mag_sq)
          {
            F -= (U_Rel*F)*U_Rel/U_Rel_mag_sq;
          }
      }

    if (dFdU)
      {
        for (unsigned int i=0; i != 3; ++i)
          for (unsigned int j=0; j != 3; ++j)
            (*dFdU)(i,j) = -(U_N(j))*U_N_unit(i);

        if (do_correction)
          {
            libmesh_not_implemented();
          }
      }

    return true;
  }

} // namespace GRINS
