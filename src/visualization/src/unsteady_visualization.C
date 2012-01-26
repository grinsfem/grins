//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "unsteady_visualization.h"

GRINS::UnsteadyVisualization::UnsteadyVisualization(const GetPot& input)
  : Visualization(input)
{
  return;
}

GRINS::UnsteadyVisualization::~UnsteadyVisualization()
{
  return;
}

void GRINS::UnsteadyVisualization::output_residual( libMesh::EquationSystems* equation_system,
						    GRINS::MultiphysicsSystem* system,
						    const unsigned int time_step )
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
  
  AutoPtr<TimeSolver> prev_time_solver(this->_system->time_solver);

  libMesh::SteadySolver* steady_solver = new libMesh::SteadySolver( *(this->_system) );

  this->_system->time_solver = AutoPtr<TimeSolver>(steady_solver);

  this->_system->assembly( true /*residual*/, false /*jacobian*/ );
  this->_system->rhs->close();

  // Swap solution with newly computed residual
  this->_system->solution->swap( *(this->_system->rhs) );
  // Update equation systems
  this->_equation_systems->update();
  
  this->dump_visualization( filename, time_step );
  
  // Now swap back and reupdate
  this->_system->solution->swap( *(this->_system->rhs) );
  this->_equation_systems->update();

  this->_system->time_solver = prev_time_solver;

  return;
}
