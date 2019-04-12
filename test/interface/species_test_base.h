//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SPECIES_TEST_BASE_H
#define GRINS_SPECIES_TEST_BASE_H

#include <vector>

// GRINS
#include "grins/physical_constants.h"

// libMesh
#include "libmesh/libmesh.h"

namespace GRINSTesting
{
  class SpeciesTestBase
  {
  public:

    SpeciesTestBase()
      : _N2_idx(libMesh::invalid_uint),
        _O2_idx(libMesh::invalid_uint),
        _O_idx(libMesh::invalid_uint),
        _N_idx(libMesh::invalid_uint),
        _NO_idx(libMesh::invalid_uint)
    {}

    libMesh::Real molar_mass( unsigned int idx )
    {
      libMesh::Real value = 0.0;

      if(idx == _N2_idx)
        value = 28.016;

      else if( idx == _O2_idx)
        value = 32.0;

      else if( idx == _NO_idx)
        value = 30.008;

      else if(idx == _O_idx)
        value = 16.0;

      else if(idx == _N_idx)
        value = 14.008;

      else
        CPPUNIT_FAIL("Invalid idx for molar_mass");

      return value;
    }

    libMesh::Real R_species( unsigned int idx )
    {
      return GRINS::Constants::R_universal/this->molar_mass(idx);
    }

  protected:

    // Species indices. Should be set by subclass at init time.
    unsigned int _N2_idx, _O2_idx, _O_idx, _N_idx, _NO_idx;

    // Species being tested
    std::vector<unsigned int> _active_species;

  };

} // end namespace GRINSTesting

#endif // GRINS_SPECIES_TEST_BASE_H
