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


#ifndef GRINS_UNSTEADY_SOLVER_H
#define GRINS_UNSTEADY_SOLVER_H

//GRINS
#include "grins/grins_solver.h"
#include "grins/adaptive_time_stepping_options.h"

//libMesh
#include "libmesh/system_norm.h"
#include "libmesh/unsteady_solver.h"

namespace GRINS
{
  class UnsteadySolver : public Solver
  {
  public:

    UnsteadySolver( const GetPot& input );
    virtual ~UnsteadySolver(){};

    virtual void solve( SolverContext& context );

  protected:

    virtual void init_time_solver(GRINS::MultiphysicsSystem* system);

    template <typename T>
    void set_theta( libMesh::UnsteadySolver* time_solver );

    //! Updates Dirichlet boundary conditions
    /*! If the Dirichlet boundary condition is nonlinear or time-dependent,
        we need to update the constraints with the new solution. */
    void update_dirichlet_bcs( SolverContext& context );

    void init_second_order_in_time_solvers( SolverContext& context );

    std::string _time_solver_name;

    unsigned int _n_timesteps;
    unsigned int _backtrack_deltat;
    double _theta;
    double _deltat;

    // Options for adaptive time solvers
    AdaptiveTimeSteppingOptions _adapt_time_step_options;

    //! Track whether is this a second order (in time) solver or not
    /*! If it is, we need to potentially initialize the acceleration */
    bool _is_second_order_in_time;

  };

  template <typename T>
  inline
  void UnsteadySolver::set_theta( libMesh::UnsteadySolver* time_solver )
  {
    T* deriv_solver = libMesh::cast_ptr<T*>(time_solver);
    deriv_solver->theta = this->_theta;
  }

} // end namespace GRINS
#endif // GRINS_UNSTEADY_SOLVER_H
