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

#ifndef GRINS_RUNNER_H
#define GRINS_RUNNER_H

// GRINS
#include "grins/simulation_initializer.h"
#include "grins/simulation_builder.h"
#include "grins/simulation.h"

// libMesh
#include "libmesh/getpot.h"

// C++
#include <iostream>

namespace GRINS
{
  //! Class to encapsulate initializing and running GRINS Simulation
  /*! This class encapsulates all of the construction, initialization, etc.
    libMesh and GRINS objects to facilitate easy construction of a GRINS-based
    program. The user only needs to construct this object, init(), and then run().
    Accessors are provided to perform auxillary functions as needed beyond a single
    Simulation run. */
  class Runner
  {
  public:

    //! Setup input objects
    /*! The constructor will only setup input parsing objects and LibMeshInit.
      The user must call init() to initialize the Simulation objects. This allows
      the user to use the underlying command line and input file objects to do steps
      that may be needed before Simulation initialization. */
    Runner( int argc, char* argv[] );

    ~Runner(){}

    const GetPot & get_input_file() const
    { return *_inputfile; }

    const GetPot & get_command_line() const
    { return _command_line; }

    Simulation & get_simulation()
    { libmesh_assert(_simulation);
      return *_simulation; }

    const Simulation & get_simulation() const
    { libmesh_assert(_simulation);
      return *_simulation; }

    const libMesh::LibMeshInit & get_libmesh_init() const
    { return *_libmesh_init; }

    //! Initialize the Simulation objects
    void init();

    //! Runs the simulation that was setup at construction time.
    void run();

  protected:

    //! Echo GRINS, libMesh version info as well as command line
    void echo_version_info( std::ostream & out, int argc, char* argv[] );

    //! Check (and error if not found) and then return GetPot input file name.
    std::string check_and_get_inputfile(int argc, char* argv[], GetPot & command_line);

    //! Check for any unused variables in GetPot input file.
    void check_for_unused_vars( const GetPot& input, bool warning_only );

    SimulationInitializer _initializer;

    SimulationBuilder _sim_builder;

    GetPot _command_line;

    std::unique_ptr<libMesh::LibMeshInit> _libmesh_init;

    std::unique_ptr<GetPot> _inputfile;

    std::unique_ptr<Simulation> _simulation;

  private:

    Runner();

  };
} // end namespace GRINS

#endif // GRINS_RUNNER_H
