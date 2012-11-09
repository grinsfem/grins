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

#include "blottner_viscosity.h"

namespace GRINS
{
  
  BlottnerViscosity::BlottnerViscosity( Real a, Real b, Real c )
    : _a(a),
      _b(b),
      _c(c)
  {
    return;
  }

  BlottnerViscosity::BlottnerViscosity( const GetPot& input )
    : _a( input("Physics/Viscosity/blottner_a", -1.0 ) ),
      _b( input("Physics/Viscosity/blottner_b", -1.0 ) ),
      _c( input("Physics/Viscosity/blottner_c", -1.0 ) )
  {
    return;
  }

  BlottnerViscosity::~BlottnerViscosity()
  {
    return;
  }

  Real BlottnerViscosity::mu( Real T )
  {
    Real logT = std::log(T);
    
    return 0.1*std::exp( (_a*logT + _b)*logT + _c );
  }

} // namespace GRINS
