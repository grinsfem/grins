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

#include "cea_curve_fit.h"

namespace GRINS
{
  CEACurveFit::CEACurveFit( const std::vector<Real>& coeffs )
    : _n_coeffs(10),
      _coefficients(coeffs)
  {
    return;
  }

  CEACurveFit::~CEACurveFit()
  {
    return;
  }

  unsigned int CEACurveFit::interval (const Real T) const
  {
    unsigned int interval = -1;

    /* CEA thermodynamic intervals are:
       [200-1,000], [1,000-6,000], [6,000-20,000] K */
    /*! \todo This could be generalized */
    if (T > 6000.)	  
      {
	interval = 2;
      }
    else if (T > 1000.)
      {
	interval =  1;
      }
    else
      {
	interval = 0;
      }

    return interval;
  }

}
