//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// GRINS
#include <memory>
#include "grins/multiphysics_sys.h"
#include "grins/mesh_builder.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/equation_systems.h"

namespace GRINSTesting
{
  //! Helper class for setting up basic GRINS::MultiphysicsSystem for unit testing
  class SystemHelper
  {
  protected:

    void setup_multiphysics_system(const std::string& filename)
    {
      _input.reset( new GetPot(filename) );
      GRINS::MeshBuilder mesh_builder;
      _mesh = mesh_builder.build( *_input, *TestCommWorld );
      _es.reset( new libMesh::EquationSystems(*_mesh) );
      _system = &_es->add_system<GRINS::MultiphysicsSystem>( "GRINS-TEST" );

      // We may not need any of these options, but this does some setup work that's needed
      // and all the options have sane defaults if they're not in the testing input file.
      _system->read_input_options( (*_input) );
    }

    void reset_all()
    {
      _input.reset();
      _es.reset(); // This will delete the system
      _mesh.reset();
    }

    std::unique_ptr<GetPot> _input;
    std::shared_ptr<libMesh::UnstructuredMesh> _mesh;
    std::unique_ptr<libMesh::EquationSystems> _es;

    // Needs to be an ordinar pointer since EquationSystems owns this
    GRINS::MultiphysicsSystem* _system;
  };

} // end namespace GRINSTesting
