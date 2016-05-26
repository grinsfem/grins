//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/homogeneous_dirichlet_bc_factory.h"

namespace GRINS
{
  // Instantiate all the AllZeroVars factories.
  HomogeneousDirichletBCFactory grins_factory_homogeneous_dirichlet("homogeneous_dirichlet");
  HomogeneousDirichletBCFactory grins_factory_no_slip("no_slip");
  HomogeneousDirichletBCFactory grins_factory_no_slip_old_style("no_slip_old_style");
  HomogeneousDirichletBCFactory grins_factory_pinned("pinned");
  HomogeneousDirichletBCFactory grins_factory_pinned_old_style("pinned_old_style");
} // end namespace GRINS
