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

#ifndef GRINS_VISUALIZATION_H
#define GRINS_VISUALIZATION_H

#include <string>
#include <vector>

// libMesh
#include "getpot.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "exodusII_io.h"
#include "vtk_io.h"

namespace GRINS
{
  class Visualization
  {
  public:

    Visualization( const GetPot& input );
    ~Visualization();

    void output();
    void output( unsigned int time_step );

    void output_residual( libMesh::EquationSystems* equation_system,
			  GRINS::MultiphysicsSystem* system );

    virtual void output_residual( libMesh::EquationSystems* equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step ) =0;

  protected:

    void dump_visualization( const std::string filename_prefix, const int time_step );

    // Visualization options
    std::string _vis_output_file_prefix;
    std::vector<std::string> _output_format;
  };
}// namespace GRINS
#endif // GRINS_VISUALIZATION_H
