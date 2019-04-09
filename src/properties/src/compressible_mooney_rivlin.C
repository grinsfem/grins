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

// This class
#include "grins/compressible_mooney_rivlin.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  CompressibleMooneyRivlin::CompressibleMooneyRivlin( const GetPot & input, const std::string & material )
    : HyperelasticStrainEnergy<CompressibleMooneyRivlin>(),
    ParameterUser("CompressibleMooneyRivlin"),
    _C1(-1),
    _C2(-1),
    _C3(-1)
  {
    this->set_parameter
      (_C1, input, "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C1", _C1);

    this->set_parameter
      (_C2, input, "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C2", _C2);

    this->set_parameter
      (_C3, input, "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C3", _C3);
  }
} // end namespace GRINS
