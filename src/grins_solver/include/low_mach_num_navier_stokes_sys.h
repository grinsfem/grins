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
// Declarations for the LowMachNumberNavierStokesSystem class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef LOW_MACH_NUM_NAVIER_STOKES_SYS_H
#define LOW_MACH_NUM_NAVIER_STOKES_SYS_H

#include <string>

#include "libmesh.h"
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
    LowMachNumberNavierStokesSystem(libMesh::EquationSystems& es,
				    const std::string& name,
				    const unsigned int number)
      : FEMSystem(es, name, number)
    {}
    
    // Destructor
    ~LowMachNumberNavierStokesSystem() {}
    
    // System initialization
    virtual void init_data ();

    // Context initialization
    virtual void init_context(libMesh::DiffContext &context);

    // Element residual and jacobian calculations
    // Time dependent parts
    virtual bool element_time_derivative (bool request_jacobian,
					  libMesh::DiffContext& context);
    
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context);
    
    // Constraint parts
    virtual bool side_constraint (bool request_jacobian,
				  libMesh::DiffContext& context);
    
    // Mass matrix part
    virtual bool mass_residual (bool request_jacobian,
				libMesh::DiffContext& context);
    
  private:

    // Indices for each variable;
    unsigned int Dummy_var;
  };

} //End namespace block

#endif // LOW_MACH_NUM_NAVIER_STOKES_SYS_H
