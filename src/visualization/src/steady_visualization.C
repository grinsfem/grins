//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/steady_visualization.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  SteadyVisualization::SteadyVisualization(const GetPot& input)
    : Visualization(input)
  {
    return;
  }

  SteadyVisualization::~SteadyVisualization()
  {
    return;
  }

  void SteadyVisualization::output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
					     MultiphysicsSystem* system,
					     const unsigned int,
					     const Real )
  {
    std::string filename = this->_vis_output_file_prefix+"_residual";

    // Idea is that this->rhs stashes the residual. Thus, when we swap
    // with the solution, we should be dumping the residual. Then, we swap
    // back once we're done outputting.

    // Swap solution with computed residual
    system->solution->swap( *(system->rhs) );
    equation_system->update();
  
    this->dump_visualization( equation_system, filename, 0.0 );
  
    // Now swap back and reupdate
    system->solution->swap( *(system->rhs) );
    equation_system->update();

    return;
  }

} // namespace GRINS
