//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of GRINS.
//
// GRINS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GRINS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GRINS.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// Declarations for the GRINSSolver class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

#include "libmesh.h"
#include "mesh.h"
#include "getpot.h"

//PB: Commented this out because I'm not sure we'll need it here.
//#include "LowMachNumberNavierStokesSystem.h"

namespace GRINS
{

  class GRINSSolver
  {
    
  public:
    GRINSSolver( const std::string application_options );
    ~GRINSSolver();

    void read_input_options( const GetPot& input );

    // get/set libMesh::Mesh
    libMesh::Mesh* get_mesh();
    void set_mesh( libMesh::Mesh* mesh );

  private:
    std::string _application_options;
    
    bool _output_vis_flag;

    libMesh::Mesh* _mesh;
  };

} //End namespace block

#endif //GRINS_SOLVER_H
