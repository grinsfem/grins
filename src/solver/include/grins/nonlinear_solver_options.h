//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_NONLINEAR_SOLVER_OPTIONS_H
#define GRINS_NONLINEAR_SOLVER_OPTIONS_H

// C++
#include <string>

// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  //! Container for nonlinear solver options
  class NonlinearSolverOptions
  {
  public:

    NonlinearSolverOptions( const GetPot & input );
    ~NonlinearSolverOptions(){};

    inline
    unsigned int max_nonlinear_iterations() const
    { return _input(_prefix+"/max_nonlinear_iterations",10); }

    inline
    libMesh::Real relative_step_tolerance() const
    { return _input(_prefix+"/relative_step_tolerance",1.0e-6); }

    inline
    libMesh::Real absolute_step_tolerance() const
    { return _input(_prefix+"/absolute_step_tolerance",0.0); }

    /* [PB]: Captured this comment from the previous location.
      _relative_residual_tolerance applies to both one of the
      stopping criteria for (nonlinear) forward solves and *the*
      stopping criterion for (linear) adjoint solves. */
    inline
    libMesh::Real relative_residual_tolerance() const
    { return _input(_prefix+"/relative_residual_tolerance",1.0e-15); }

    inline
    libMesh::Real absolute_residual_tolerance() const
    { return _input(_prefix+"/absolute_residual_tolerance",0.0); }

    inline
    bool continue_after_backtrack_failure() const
    { return _input(_prefix+"/continue_after_backtrack_failure",false); }

    inline
    bool continue_after_max_iterations() const
    { return _input(_prefix+"/continue_after_max_iterations",false); }

    inline
    bool require_residual_reduction() const
    { return _input(_prefix+"/require_residual_reduction",true); }

  protected:

    const GetPot & _input;

    std::string _prefix;

  private:

    NonlinearSolverOptions();

  };

} // end namespace GRINS

#endif // GRINS_NONLINEAR_SOLVER_OPTIONS_H
