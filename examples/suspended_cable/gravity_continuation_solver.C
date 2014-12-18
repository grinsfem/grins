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
#include "gravity_continuation_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"
#include "grins/elastic_cable_constant_gravity.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  GravityContinuationSolver::GravityContinuationSolver( const GetPot& input )
    : SteadySolver(input)
  {
    if( !input.have_variable("SolverOptions/GravityContinuation/final_gravity") )
      {
        std::cerr << "Error: Did not find final_gravity value for GravityContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/GravityContinuation/final_gravity" << std::endl;
        libmesh_error();
      }
    libMesh::Real final_gravity_value = input("SolverOptions/GravityContinuation/final_gravity", 0.0);

    libMesh::Real initial_gravity_value = input("SolverOptions/GravityContinuation/initial_gravity", 0.0);



    if( !input.have_variable("SolverOptions/GravityContinuation/n_increments") )
      {
        std::cerr << "Error: Did not find n_increments value for GravityContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/GravityContinuation/n_increments" << std::endl;
        libmesh_error();
      }
    unsigned int n_increments = input("SolverOptions/GravityContinuation/n_increments", 1);

    _gravity_values.resize(n_increments);

    libMesh::Real increment = (final_gravity_value - initial_gravity_value)/n_increments;

    for( unsigned int i = 0; i < n_increments; i++ )
      {
    	_gravity_values[i] = (i+1)*increment + initial_gravity_value;
      }

    return;
  }

  GravityContinuationSolver::~GravityContinuationSolver()
  {
    return;
  }

  void GravityContinuationSolver::solve( SolverContext& context )
  {
    for( unsigned int s = 0; s < _gravity_values.size(); s++ )
      {

        libMesh::Real gravity = _gravity_values[s];

        std::cout << "==========================================================" << std::endl
                  << "   Gravity step " << s  << ", Gravity = " << gravity << std::endl
                  << "==========================================================" << std::endl;

        this->increment_gravity(*(context.system), gravity);

        context.system->solve();

        if( context.output_vis )
          {
            context.postprocessing->update_quantities( *(context.equation_system) );
            context.vis->output( context.equation_system, s, gravity );
          }
      }

    return;
  }

  void GravityContinuationSolver::increment_gravity( GRINS::MultiphysicsSystem& system,
                                                       libMesh::Real gravity )
  {
    // Get Physics class and cast
    std::tr1::shared_ptr<GRINS::Physics> raw_physics = system.get_physics(elastic_cable_constant_gravity);
    ElasticCableConstantGravity& physics = libMesh::cast_ref<ElasticCableConstantGravity&>( *(raw_physics.get()) );

    physics.reset_gravity(gravity);

    return;
  }

} // end namespace GRINS
