//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

#include "config.h"

#include "getpot.h"
#include "libmesh.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "diff_solver.h"
#include "equation_systems.h"
#include "euler_solver.h"
#include "steady_solver.h"

#ifdef HAVE_GRVY
#include "grvy.h" // GRVY timers
#endif

namespace GRINS
{
  // We template around the system type
  template< class T >
  class Solver
  {
    
  public:
    Solver();
    ~Solver();

    void read_input_options( const GetPot& input );

    // get/set pair for mesh object
    libMesh::Mesh* get_mesh();
    void set_mesh( libMesh::Mesh* mesh );

    //TODO: Should we have these return error codes?
    void initialize_system( std::string system_name, GetPot& input );
    void solve( );

    void output_visualization();
    void output_visualization( unsigned int time_step );

    T* get_system();
    
#ifdef USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  private:
    
    // Linear/Nonlinear solver options
    unsigned int _max_nonlinear_iterations;
    double _relative_step_tolerance;
    double _absolute_step_tolerance;
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
    std::vector<std::string> _output_format;

    // Mesh/Solver related objects
    libMesh::Mesh* _mesh;
    libMesh::EquationSystems* _equation_systems;
    T* _system;

    // Variables to track state of system
    bool _system_initialized;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _solver_quiet;
    bool _solver_verbose;
    bool _print_equation_system_info;

    void dump_visualization( std::string filename_prefix );
    void set_solver_options( libMesh::DiffSolver& solver );

#ifdef USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  };

} //End namespace block

#endif //GRINS_SOLVER_H
