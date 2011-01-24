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

namespace GRINS
{

  // FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
  // but we must specify element residuals
  class LowMachNumberNavierStokesSystem : public FEMSystem
  {
    
  public:
    // Constructor
    LowMachNumberNavierStokesSystem(EquationSystems& es,
				    const std::string& name,
				    const unsigned int number)
      : FEMSystem(es, name, number)
    {}
    
    // Destructor
    ~LowMachNumberNavierStokesSystem() {}
    
    void set_application( const std::string application_options );
    
  private:
    std::string _application_options;
  };
  
} //End namespace block

#endif // LOW_MACH_NUM_NAVIER_STOKES_SYS_H
