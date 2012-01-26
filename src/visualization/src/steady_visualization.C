//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "steady_visualization.h"

GRINS::SteadyVisualization::SteadyVisualization(const GetPot& input)
  : Visualization(input)
{
  return;
}

GRINS::SteadyVisualization::~SteadyVisualization()
{
  return;
}

void GRINS::SteadyVisualization::output_residual( libMesh::EquationSystems* equation_system,
						  GRINS::MultiphysicsSystem* system,
						  const unsigned int time_step )
{
  std::stringstream suffix;
  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix+"_residual";

  filename+="."+suffix.str();

  // Idea is that this->rhs stashes the residual. Thus, when we swap
  // with the solution, we should be dumping the residual. Then, we swap
  // back once we're done outputting.

  // Swap solution with computed residual
  system->solution->swap( *(system->rhs) );
  equation_system->update();
  
  this->dump_visualization( filename, time_step );
  
  // Now swap back and reupdate
  system->solution->swap( *(system->rhs) );
  equation_system->update();

  return;
}
