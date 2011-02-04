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
// Declarations for the GRINS::Solver class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

#include "getpot.h"
#include "libmesh.h"
#include "mesh.h"
#include "diff_solver.h"
#include "equation_systems.h"
#include "euler_solver.h"
#include "steady_solver.h"

namespace GRINS
{
  // We template around the system type
  template< class T >
  class Solver
  {
    
  public:
    Solver( const std::string application_options );
    ~Solver();

    void read_input_options( const GetPot& input );

    // get/set pair for mesh object
    libMesh::Mesh* get_mesh();
    void set_mesh( libMesh::Mesh* mesh );

    //TODO: Should we have these return error codes?
    void initialize_system();
    void solve( );

    void output_visualization();
    void output_visualization( unsigned int time_step );

  private:
    std::string _application_options;
    
    // Linear/Nonlinear solver options
    bool _solver_quiet;
    unsigned int _max_nonlinear_iterations;
    double _relative_step_tolerance;
    double _relative_residual_tolerance;
    double _absolute_residual_tolerance;
    unsigned int _max_linear_iterations;
    double _initial_linear_tolerance;

    // Unsteady solver options
    bool _transient;
    double _theta;
    unsigned int _n_timesteps;
    double _deltat;

    // Visualization options
    bool _output_vis_time_series;
    std::string _vis_output_file_prefix;
    std::string _output_format;

    // Mesh/Solver related objects
    libMesh::Mesh* _mesh;
    libMesh::EquationSystems* _equation_systems;
    T* _system;

    // Variables to track state of system
    bool _system_initialized;

    void dump_visualization( std::string filename_prefix );
    void set_solver_options( libMesh::DiffSolver& solver );
  };

} //End namespace block

#endif //GRINS_SOLVER_H
