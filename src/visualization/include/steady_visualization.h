//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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

#ifndef GRINS_STEADY_VISUALIZATION_H
#define GRINS_STEADY_VISUALIZATION_H

#include "visualization.h"

namespace GRINS
{
  class SteadyVisualization : public Visualization
  {
  public:

    SteadyVisualization(const GetPot& input);
    ~SteadyVisualization();

    virtual void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step,
				  const Real time );

  protected:

  };
}
#endif
