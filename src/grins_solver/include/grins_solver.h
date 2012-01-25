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

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

// GRINS
#include "config.h"
#include "multiphysics_sys.h"

// libMesh
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
  class Solver
  {
    
  public:
    Solver();
    ~Solver();

    virtual void read_input_options( const GetPot& input );

    void initialize( GetPot& input, 
		     libMesh::EquationSystems equation_system,
		     GRINS::PhysicsList& physics_list );
    
    virtual void solve()=0;
    virtual void init_time_solver()=0;

    void output_visualization();
    void output_visualization( unsigned int time_step );

    void output_residual_vis( const unsigned int time_step );
    void output_unsteady_residual_vis( const unsigned int time_step );
    
#ifdef USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:
    
    // Linear/Nonlinear solver options
    unsigned int _max_nonlinear_iterations;
    double _relative_step_tolerance;
    double _absolute_step_tolerance;
    double _relative_residual_tolerance;
    double _absolute_residual_tolerance;
    unsigned int _max_linear_iterations;
    double _initial_linear_tolerance;

    // Visualization options
    bool _output_vis_time_series;
    std::string _vis_output_file_prefix;
    std::vector<std::string> _output_format;
    bool _output_residual;
    bool _output_unsteady_residual;

    GRINS::MultiphysicsSystem* _system;

    // Screen display options
    bool _solver_quiet;
    bool _solver_verbose;
    
    void dump_visualization( const std::string filename_prefix, const int time_step );
    void set_solver_options( libMesh::DiffSolver& solver );

#ifdef USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  };

} //End namespace block

#endif //GRINS_SOLVER_H
