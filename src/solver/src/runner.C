//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/runner.h"

// libMesh
#include "libmesh/parallel.h"

// C++
#include <iomanip>
#include <sstream>

namespace GRINS
{
  Runner::Runner( int argc, char* argv[] )
    : _initializer(),
      _sim_builder(),
      _command_line(argc,argv),
      _libmesh_init( new libMesh::LibMeshInit(argc,argv) )
  {
    this->echo_version_info( libMesh::out, argc, argv );

    // Grab inputfile name and setup GetPot input file
    std::string inputfile_name = this->check_and_get_inputfile(argc,argv);
    _inputfile.reset( new GetPot(inputfile_name) );

    // Allow command line options to override the inputfile
    _inputfile->parse_command_line(argc,argv);
  }

  void Runner::init()
  {
    // Initialize Simulation
    _simulation.reset( new Simulation( *_inputfile,
                                       _command_line,
                                       _sim_builder,
                                       _libmesh_init->comm() ) );
  }

  void Runner::echo_version_info( std::ostream & out, int argc, char* argv[] )
  {
    out << "==========================================================" << std::endl;
    out << "GRINS Version: " << GRINS_BUILD_VERSION << std::endl
        << "libMesh Version: " << LIBMESH_BUILD_VERSION << std::endl
        << "Running with command:\n";

    for (int i=0; i != argc; ++i)
      out << argv[i] << ' ';

    out << std::endl
        << "==========================================================" << std::endl;
  }

  std::string Runner::check_and_get_inputfile(int argc, char* argv[])
  {
    if( argc < 2 )
    {
      std::stringstream error_msg;
      error_msg << "ERROR: Found only 1 command line argument, but was expecting an inputfile name!"
                << std::endl
                << "       Please specify the name of the inputfile on the command line." << std::endl;
      libmesh_error_msg(error_msg.str());
    }

    std::string inputfile_name = argv[1];

    std::ifstream i(inputfile_name.c_str());
    if (!i)
      {
        std::string error_msg = "Error: Could not read from input file "+inputfile_name+"!\n";
        libmesh_error_msg(error_msg);
      }

    return inputfile_name;
  }

} // end namespace GRINS
