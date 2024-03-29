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

#ifndef GRINS_DISPLACEMENT_CONTINUATION_SOLVER_H
#define GRINS_DISPLACEMENT_CONTINUATION_SOLVER_H

//GRINS
#include "grins/steady_solver.h"

namespace GRINS
{
  class DisplacementContinuationSolver : public SteadySolver
  {
  public:

    DisplacementContinuationSolver( const GetPot& input );

    virtual ~DisplacementContinuationSolver() = default;

    virtual void initialize( const GetPot& input,
                             std::shared_ptr<libMesh::EquationSystems> equation_system,
                             GRINS::MultiphysicsSystem* system ) override;

    virtual void solve( SolverContext& context ) override;

  protected:

    void increment_displacement( GRINS::MultiphysicsSystem& system,
                                 libMesh::EquationSystems& equation_system,
                                 const libMesh::Real displacement );

    //! Boundary on which we want to increment the displacement
    libMesh::boundary_id_type _bc_id;

    //! Cache index into libMesh::DirichletBoundaries
    /*!
     * libMesh::DirichletBoundaries subclasses std::vector. So, we
     * cache the index to the particular libMesh::DirichletBoundary
     * we want. This way, we search at the beginning and reuse.
     */
    unsigned int _bc_index;

    std::vector<libMesh::Real> _displacements;

  };
} // namespace GRINS
#endif // GRINS_DISPLACEMENT_CONTINUATION_SOLVER_H
