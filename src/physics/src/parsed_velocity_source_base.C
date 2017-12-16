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
#include "grins/parsed_velocity_source_base.h"

// GRINS
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{

  template<class Mu>
  ParsedVelocitySourceBase<Mu>::ParsedVelocitySourceBase( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    _input(input)
  {
  }

  template<class Mu>
  ParsedVelocitySourceBase<Mu>::~ParsedVelocitySourceBase()
  {
    return;
  }

  template<class Mu>
  void ParsedVelocitySourceBase<Mu>::set_time_evolving_vars ( libMesh::FEMSystem* system)
  {
    std::string base_physics_name = "ParsedVelocitySource";

    libMesh::ParsedFEMFunction<libMesh::Number> *vsf
      (new libMesh::ParsedFEMFunction<libMesh::Number> (*system, ""));
    this->velocity_source_function.reset(vsf);

    this->set_parameter(*vsf, _input,
                        "Physics/"+base_physics_name+"/source_function",
                        "DIE!");
  }

  template<class Mu>
  bool ParsedVelocitySourceBase<Mu>::compute_force
  ( const libMesh::Point& point,
    const libMesh::Real time,
    const AssemblyContext &c,
    libMesh::NumberVectorValue& F,
    libMesh::NumberTensorValue *dFdU)
  {
    libmesh_assert(velocity_source_function.get());

    libMesh::DenseVector<libMesh::Number> output_vec(3);

    (*velocity_source_function)(c, point, time,
                                output_vec);

    F(0) = output_vec(0);
    F(1) = output_vec(1);
    F(2) = output_vec(2);

    if (dFdU)
      for (unsigned int i=0; i != 3; ++i)
        for (unsigned int j=0; j != 3; ++j)
          (*dFdU)(i,j) = 0; // FIXME

    if (F(0) || F(1) || F(2))
      return true;

    return false;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(ParsedVelocitySourceBase);
