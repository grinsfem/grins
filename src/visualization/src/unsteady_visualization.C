//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "unsteady_visualization.h"

namespace GRINS
{

  UnsteadyVisualization::UnsteadyVisualization(const GetPot& input)
    : Visualization(input)
  {
    return;
  }

  UnsteadyVisualization::~UnsteadyVisualization()
  {
    return;
  }

  void UnsteadyVisualization::output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
					       MultiphysicsSystem* system,
					       const unsigned int time_step,
					       const Real time )
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
  
    libMesh::AutoPtr<TimeSolver> prev_time_solver(system->time_solver);

    libMesh::SteadySolver* steady_solver = new libMesh::SteadySolver( *(system) );

    system->time_solver = AutoPtr<TimeSolver>(steady_solver);

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

    system->time_solver = prev_time_solver;

    return;
  }

} // namespace GRINS
