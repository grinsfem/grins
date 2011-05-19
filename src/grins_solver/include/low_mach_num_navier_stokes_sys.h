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
// Declarations for the LowMachNumberNavierStokesSystem class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef LOW_MACH_NUM_NAVIER_STOKES_SYS_H
#define LOW_MACH_NUM_NAVIER_STOKES_SYS_H

#include <string>

#include "libmesh.h"
#include "getpot.h"
#include "fem_system.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"

// DiffSystem framework files
#include "fem_system.h"
#include "fem_context.h"

namespace GRINS
{

  // FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
  // but we must specify element residuals
  class LowMachNumberNavierStokesSystem : public FEMSystem
  {
    
  public:
    // Constructor
    LowMachNumberNavierStokesSystem( libMesh::EquationSystems& es,
				     const std::string& name,
				     const unsigned int number )
      : FEMSystem(es, name, number)
    {}
    
    // Destructor
    ~LowMachNumberNavierStokesSystem() {}
    
    virtual void read_input_options( GetPot& input );

    // System initialization
    virtual void init_data();

    // Context initialization
    virtual void init_context( libMesh::DiffContext &context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context );
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context );
    
    // Constraint part(s)
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context );
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context );
    
    // Mass matrix part(s)
    /* virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context ); */
    
  private:

    // Indices for each variable;
    unsigned int u_var;
    unsigned int v_var;
    unsigned int w_var;
    unsigned int p_var;
    unsigned int T_var;

    // Returns the value of a forcing function at point pt_xyz.
    // This value depends on which option is set.
    // TODO: any other option to return? other than libMesh::Point?
    libMesh::Point forcing(const libMesh::Point& pt_xyz);
  };

} //End namespace block

#endif // LOW_MACH_NUM_NAVIER_STOKES_SYS_H
