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

#ifndef GRINS_ANTIOCH_TEST_BASE_H
#define GRINS_ANTIOCH_TEST_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/getpot.h"

namespace GRINSTesting
{
  class AntiochTestBase
  {
  public:

    void init_antioch(const std::string& input_file, const std::string& material_name)
    {
      GetPot input(input_file);

      _antioch_mixture.reset( new GRINS::AntiochMixture(input,material_name) );
    }

  protected:

    libMesh::UniquePtr<GRINS::AntiochMixture> _antioch_mixture;
  };

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_TEST_BASE_H
