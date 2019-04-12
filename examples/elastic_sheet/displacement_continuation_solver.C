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

// This class
#include "displacement_continuation_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/const_function.h"
#include "libmesh/composite_function.h"

namespace GRINS
{

  DisplacementContinuationSolver::DisplacementContinuationSolver( const GetPot& input )
    : SteadySolver(input),
      _bc_id( input("SolverOptions/DisplacementContinuation/boundary", -1) )
  {
    if( !input.have_variable("SolverOptions/DisplacementContinuation/boundary") )
      {
        std::cerr << "Error: Did not find boundary id for DisplacementContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/DisplacementContinuation/boundary in input." << std::endl;
        libmesh_error();
      }

    if( !input.have_variable("SolverOptions/DisplacementContinuation/final_displacement") )
      {
        std::cerr << "Error: Did not find final_displacement value for DisplacementContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/DisplacementContinuation/final_displacement" << std::endl;
        libmesh_error();
      }
    libMesh::Real disp_value = input("SolverOptions/DisplacementContinuation/final_displacement", 0.0);


    if( !input.have_variable("SolverOptions/DisplacementContinuation/n_increments") )
      {
        std::cerr << "Error: Did not find final_displacement value for DisplacementContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/DisplacementContinuation/n_increments" << std::endl;
        libmesh_error();
      }
    unsigned int n_increments = input("SolverOptions/DisplacementContinuation/n_increments", 1);

    _displacements.resize(n_increments);

    libMesh::Real increment = disp_value/n_increments;

    for( unsigned int i = 0; i < n_increments; i++ )
      {
        _displacements[i] = (i+1)*increment;
      }

    return;
  }

  DisplacementContinuationSolver::~DisplacementContinuationSolver()
  {
    return;
  }

  void DisplacementContinuationSolver::initialize( const GetPot& input,
                                                   std::shared_ptr<libMesh::EquationSystems> equation_system,
                                                   GRINS::MultiphysicsSystem* system )
  {
    // First init everything on the base class side, which will reinit equation_system
    Solver::initialize(input,equation_system,system);

    // Now search for the index for the DirichletBoundary we want to increment
    libMesh::DirichletBoundaries* d_vector = system->get_dof_map().get_dirichlet_boundaries();

    bool found_bc_id = false;
    for( libMesh::DirichletBoundaries::const_iterator it = d_vector->begin(); it != d_vector->end(); ++it )
      {
        if( (*it)->b.find(_bc_id) != (*it)->b.end() )
          {
            found_bc_id = true;
            _bc_index = it - d_vector->begin();
            // We're assuming that there's only one boundary
            break;
          }
      }

    if( !found_bc_id )
      {
        std::cerr << "Error: Did not find prescribed boundary for DisplacementContinuationSolver" << std::endl
                  << "       Was searching for bc_id = " << _bc_id << std::endl;
        libmesh_error();
      }

    return;
  }

  void DisplacementContinuationSolver::solve( SolverContext& context )
  {
    for( unsigned int s = 0; s < _displacements.size(); s++ )
      {

        libMesh::Real disp = _displacements[s];

        std::cout << "==========================================================" << std::endl
                  << "   Displacement step " << s  << ", displacement = " << disp << std::endl
                  << "==========================================================" << std::endl;

        this->increment_displacement(*(context.system), *(context.equation_system), disp);

        context.system->solve();

        if( context.output_vis )
          {
            context.postprocessing->update_quantities( *(context.equation_system) );
            context.vis->output( context.equation_system, s, disp );
          }
      }

    return;
  }

  void DisplacementContinuationSolver::increment_displacement( GRINS::MultiphysicsSystem& system,
                                                               libMesh::EquationSystems& equation_system,
                                                               const libMesh::Real displacement )
  {
    // Get DirichetBoundaries vector.
    libMesh::DirichletBoundaries* d_vector = system.get_dof_map().get_dirichlet_boundaries();

    // Get the DirichletBoundary we want
    libMesh::DirichletBoundary* dirichlet = (*d_vector)[_bc_index];

    // Kill the old FunctionBase object and put in our new one.
    libMesh::FunctionBase<libMesh::Real>* composite_func_ptr = new libMesh::CompositeFunction<libMesh::Real>;
    libMesh::CompositeFunction<libMesh::Real>& composite_func = libMesh::cast_ref<libMesh::CompositeFunction<libMesh::Real>&>( *composite_func_ptr );

    std::vector<VariableIndex> var_idx(1,0); // Hardcoding to Ux displacement component
    composite_func.attach_subfunction( libMesh::ConstFunction<libMesh::Real>(displacement), var_idx );

    // DirichletBoundary now takes ownership of the pointer
    dirichlet->f.reset(composite_func_ptr);

    // Need to reinit system
    equation_system.reinit();

    return;
  }


} // end namespace GRINS
