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
#include "grins/unsteady_visualization.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/steady_solver.h"

namespace GRINS
{

  UnsteadyVisualization::UnsteadyVisualization
  ( const GetPot& input,
    const libMesh::Parallel::Communicator &comm )
    : Visualization(input, comm)
  {
    return;
  }

  UnsteadyVisualization::~UnsteadyVisualization()
  {
    return;
  }

  void UnsteadyVisualization::output_residual
  ( std::shared_ptr<libMesh::EquationSystems> equation_system,
    MultiphysicsSystem* system,
    const unsigned int time_step,
    const libMesh::Real time )
  {
    std::stringstream suffix;
    suffix << time_step;

    std::string filename = this->_vis_output_file_prefix+"_unsteady_residual";

    filename+="."+suffix.str();

    // For the unsteady residual, we just want to evaluate F(u) from
    // dU/dt = F(u). What we do is swap out the time solver to a
    // SteadySolver and reassemble the residual. Then, we'll need to swap
    // the solution and the rhs vector stashed in the system. Once we're done,
    // we'll reset the time solver pointer back to the original guy.

    libMesh::TimeSolver* prev_time_solver = system->time_solver.get();

    libMesh::SteadySolver* steady_solver = new libMesh::SteadySolver( *(system) );

    system->time_solver = std::unique_ptr<libMesh::TimeSolver>(steady_solver);

    system->assembly( true /*residual*/, false /*jacobian*/ );
    system->rhs->close();

    // Swap solution with newly computed residual
    system->solution->swap( *(system->rhs) );
    // Update equation systems
    equation_system->update();

    this->dump_visualization( equation_system, filename, time );

    // Now swap back and reupdate
    system->solution->swap( *(system->rhs) );
    equation_system->update();

    system->time_solver.reset(prev_time_solver);
  }

  void UnsteadyVisualization::output_residual_sensitivities
  (std::shared_ptr<libMesh::EquationSystems> equation_system,
   MultiphysicsSystem* system,
   const libMesh::ParameterVector & params,
   const unsigned int time_step,
   const libMesh::Real time )
  {
    for (unsigned int p=0; p != params.size(); ++p)
      {
        std::stringstream suffix;
        suffix << time_step;

        std::stringstream pstr;
        pstr << p;

        std::string filename =
          this->_vis_output_file_prefix + "_unsteady_dRdp" +
          pstr.str() + '.' + suffix.str();

        // Swap solution with precomputed sensitivity rhs
        system->solution->swap(system->get_sensitivity_rhs(p));
        equation_system->update();

        this->dump_visualization( equation_system, filename, time );

        // Now swap back and reupdate
        system->solution->swap(system->get_sensitivity_rhs(p));
        equation_system->update();
      }
  }

  void UnsteadyVisualization::output_adjoint
  ( std::shared_ptr<libMesh::EquationSystems> equation_system,
    MultiphysicsSystem* system,
    const unsigned int time_step,
    const libMesh::Real time )
  {
    std::stringstream suffix;
    suffix << time_step;

    const libMesh::DifferentiableQoI* raw_qoi = system->get_qoi();
    const CompositeQoI* qoi = dynamic_cast<const CompositeQoI*>( raw_qoi );

    unsigned int n_qois = qoi->n_qois();

    for( unsigned int q = 0; q < n_qois; q++ )
      {
        libMesh::NumericVector<libMesh::Number>& dual_solution = system->get_adjoint_solution(q);

        const std::string& qoi_name = qoi->get_qoi(q).name();
        std::string filename = this->_vis_output_file_prefix+"_unsteady_adjoint_"+qoi_name;
        filename+="."+suffix.str();

        system->solution->swap( dual_solution );
        equation_system->update();

        this->dump_visualization( equation_system, filename, time );

        // Now swap back and reupdate
        system->solution->swap( dual_solution );
        equation_system->update();
      }
  }

  void UnsteadyVisualization::output_solution_sensitivities
  (std::shared_ptr<libMesh::EquationSystems> equation_system,
   MultiphysicsSystem* system,
   const libMesh::ParameterVector & params,
   const unsigned int time_step,
   const libMesh::Real time )
  {
    for (unsigned int p=0; p != params.size(); ++p)
      {
        std::stringstream suffix;
        suffix << time_step;

        std::stringstream pstr;
        pstr << p;

        std::string filename =
          this->_vis_output_file_prefix + "_unsteady_dudp" +
          pstr.str() + '.' + suffix.str();

        // Swap solution with precomputed sensitivity solution
        system->solution->swap(system->get_sensitivity_solution(p));
        equation_system->update();

        this->dump_visualization( equation_system, filename, time );

        // Now swap back and reupdate
        system->solution->swap(system->get_sensitivity_solution(p));
        equation_system->update();
      }
  }


} // namespace GRINS
