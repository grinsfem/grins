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

#ifndef GRINS_VISUALIZATION_H
#define GRINS_VISUALIZATION_H

#include <string>
#include <vector>
#include "boost/tr1/memory.hpp"

// libMesh
#include "numeric_vector.h"
#include "equation_systems.h"
#include "getpot.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "exodusII_io.h"
#include "vtk_io.h"
#include "steady_solver.h"

// GRINS
#include "multiphysics_sys.h"

namespace GRINS
{
  class Visualization
  {
  public:

    Visualization( const GetPot& input );
    ~Visualization();

    void output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system );
    void output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
		 const unsigned int time_step, const Real time );

    void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			  GRINS::MultiphysicsSystem* system );

    virtual void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step, const Real time ) =0;

    void dump_visualization( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			     const std::string& filename_prefix, const Real time );
    
  protected:

    // Visualization options
    std::string _vis_output_file_prefix;
    std::vector<std::string> _output_format;
  };
}// namespace GRINS
#endif // GRINS_VISUALIZATION_H
