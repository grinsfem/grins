//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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

#ifndef GRINS_REACTING_FLOW_CACHE_H
#define GRINS_REACTING_FLOW_CACHE_H

// C++
#include <vector>

// libMesh
#include "libmesh_common.h"

namespace GRINS
{
  class ReactingFlowCache
  {
  public:

    ReactingFlowCache( Real T, Real P, std::vector<Real>& mass_fractions );
    ~ReactingFlowCache();

    Real T() const
    { return _T; };

    Real P() const
    { return _P; }

    const std::vector<Real>& mass_fractions() const
    { return _mass_fractions; }
    
  protected:

    const Real _T;
    const Real _P;
    const std::vector<Real> _mass_fractions;

  private:
    
    ReactingFlowCache();

  };

} // namespace GRINS

#endif //GRINS_REACTING_FLOW_CACHE_H
