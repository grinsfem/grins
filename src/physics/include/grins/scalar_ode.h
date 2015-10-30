//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_SCALAR_ODE_H
#define GRINS_SCALAR_ODE_H

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/cached_values.h"
#include "grins/inc_navier_stokes_base.h"

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/getpot.h"

// C++
#include <string>

namespace GRINS
{

  //! Physics class for arbitrary scalar-valued ODEs
  /*
    This physics class allows ODEs specified in ParsedFEMFunction
    config file arguments to be solved.
   */
  class ScalarODE : public Physics
  {
  public:

    ScalarODE( const std::string& physics_name, const GetPot& input );

    ~ScalarODE();

    //! Initialization of variables
    /*!
      Add scalar variable(s) to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets scalar variable(s) to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Prepare the context for evaluations
    virtual void init_context( AssemblyContext& context );
    
    // residual and jacobian calculations

    // User-specified ODE(s)
    virtual void nonlocal_time_derivative ( bool compute_jacobian,
				            AssemblyContext& context,
				            CachedValues& cache );

    // User-specified constraint equation
    virtual void nonlocal_constraint ( bool compute_jacobian,
				       AssemblyContext& context,
				       CachedValues& cache );

    // User-specified (or default "s'") mass term
    virtual void nonlocal_mass_residual ( bool compute_jacobian,
				          AssemblyContext& context,
				          CachedValues& cache );

    VariableIndex scalar_ode_var() const { return _scalar_ode_var; }

  private:

    // ParsedFEMFunctions evaluating the mass, time derivative, and
    // constraint components of an ODE.
    libMesh::AutoPtr<libMesh::FEMFunctionBase<libMesh::Number> >
      time_deriv_function,
      constraint_function,
      mass_residual_function;

    // Number of components of the scalar solution variable; defaults
    // to 1.
    libMesh::Number _order;

    // Perturbation to use for finite differencing of functions
    libMesh::Number _epsilon;

    VariableIndex _scalar_ode_var; /* Index for turbine speed scalar */

    std::string _scalar_ode_var_name;

    const GetPot & _input;

    ScalarODE();
  };

} // end namespace block

#endif // GRINS_SCALAR_ODE_H
