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
#include "pressure_continuation_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"
#include "grins/elastic_membrane_constant_pressure.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  PressureContinuationSolver::PressureContinuationSolver( const GetPot& input )
    : SteadySolver(input)
  {
    if( !input.have_variable("SolverOptions/PressureContinuation/final_pressure") )
      {
        std::cerr << "Error: Did not find final_pressure value for PressureContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/PressureContinuation/final_pressure" << std::endl;
        libmesh_error();
      }
    libMesh::Real final_pressure_value = input("SolverOptions/PressureContinuation/final_pressure", 0.0);

    libMesh::Real initial_pressure_value = input("SolverOptions/PressureContinuation/initial_pressure", 0.0);



    if( !input.have_variable("SolverOptions/PressureContinuation/n_increments") )
      {
        std::cerr << "Error: Did not find n_increments value for PressureContinuationSolver" << std::endl
                  << "       Must specify SolverOptions/PressureContinuation/n_increments" << std::endl;
        libmesh_error();
      }
    unsigned int n_increments = input("SolverOptions/PressureContinuation/n_increments", 1);

    _pressure_values.resize(n_increments);

    libMesh::Real increment = (final_pressure_value - initial_pressure_value)/n_increments;

    for( unsigned int i = 0; i < n_increments; i++ )
      {
        _pressure_values[i] = (i+1)*increment + initial_pressure_value;
      }

    return;
  }

  PressureContinuationSolver::~PressureContinuationSolver()
  {
    return;
  }

  void PressureContinuationSolver::solve( SolverContext& context )
  {
    for( unsigned int s = 0; s < _pressure_values.size(); s++ )
      {

        libMesh::Real pressure = _pressure_values[s];

        std::cout << "==========================================================" << std::endl
                  << "   Pressure step " << s  << ", pressure = " << pressure << std::endl
                  << "==========================================================" << std::endl;

        this->increment_pressure(*(context.system), pressure);

        context.system->solve();

        if( context.output_vis )
          {
            context.postprocessing->update_quantities( *(context.equation_system) );
            context.vis->output( context.equation_system, s, pressure );
          }
      }

    return;
  }

  void PressureContinuationSolver::increment_pressure( GRINS::MultiphysicsSystem& system,
                                                       libMesh::Real pressure )
  {
    // Get Physics class and cast
    std::shared_ptr<GRINS::Physics> raw_physics = system.get_physics(GRINS::PhysicsNaming::elastic_membrane_constant_pressure());
    ElasticMembraneConstantPressure& physics = libMesh::cast_ref<ElasticMembraneConstantPressure&>( *(raw_physics.get()) );

    physics.reset_pressure(pressure);

    return;
  }

} // end namespace GRINS
