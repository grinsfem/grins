//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins/arrhenius_rate.h"

namespace GRINS
{
  ArrheniusRate::ArrheniusRate(const libMesh::Real Cf, const libMesh::Real eta, const libMesh::Real Ea)
    : _Cf(Cf),
      _eta(eta),
      _Ea(Ea)
  {
    return;
  }

  ArrheniusRate::~ArrheniusRate()
  {
    return;
  }

  void ArrheniusRate::print(std::ostream& os) const
  {
    os << _Cf;
    if (_eta != 0.) os << "*T^" << _eta;
    os << "*exp(-" << _Ea << "/T)";

    return;
  }

  std::ostream& operator<<(std::ostream& os, const ArrheniusRate& rate)
  {
    rate.print(os);
    return os;
  }

} // namespace GRINS
