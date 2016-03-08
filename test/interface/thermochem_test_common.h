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

#ifndef GRINS_THERMOCHEM_TEST_COMMON_H
#define GRINS_THERMOCHEM_TEST_COMMON_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>

#include <vector>

#include "libmesh/libmesh_common.h"

namespace GRINSTesting
{
  class ThermochemTestCommon
  {
  public:

    static libMesh::Real compute_mass_frac_mixture_prop( const std::vector<libMesh::Real>& properties,
                                                         const std::vector<libMesh::Real>& mass_fracs )
    {
      CPPUNIT_ASSERT_EQUAL( properties.size(), mass_fracs.size() );

      unsigned int size = properties.size();

      libMesh::Real mixed_value = 0.0;
      for( unsigned int s = 0; s < size; s++ )
        mixed_value += mass_fracs[s]*properties[s];

      return mixed_value;
    }

    static libMesh::Real arrhenius_rate( libMesh::Real A, /* preexponential */
                                         libMesh::Real b, /* temp. exponent */
                                         libMesh::Real Ea, /* Activation energy */
                                         libMesh::Real T )
    {
      libMesh::Real RT = GRINS::Constants::R_universal*T;

      return A*std::pow(T,b)*std::exp(-Ea/RT);
    }

  };

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif // GRINS_THERMOCHEM_TEST_COMMON_H
