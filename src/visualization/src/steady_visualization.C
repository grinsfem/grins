//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
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
