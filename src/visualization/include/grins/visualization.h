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

#ifndef GRINS_VISUALIZATION_H
#define GRINS_VISUALIZATION_H

// C++
#include <string>
#include <vector>
#include "boost/tr1/memory.hpp"

// libMesh
#include "libmesh/equation_systems.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  class Visualization
  {
  public:

    Visualization( const GetPot& input );
    ~Visualization();

    void output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system );
    void output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
		 const unsigned int time_step, const libMesh::Real time );

    void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			  GRINS::MultiphysicsSystem* system );

    virtual void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step, const libMesh::Real time ) =0;

    void dump_visualization( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			     const std::string& filename_prefix, const libMesh::Real time );
    
  protected:

    // Visualization options
    std::string _vis_output_file_prefix;
    std::vector<std::string> _output_format;
  };
}// namespace GRINS
#endif // GRINS_VISUALIZATION_H
