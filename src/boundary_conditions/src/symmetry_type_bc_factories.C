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
#include "grins/symmetry_type_bc_factories.h"

namespace GRINS
{
  // Instantiate all the SymmetryType factories.
  YZSymmetryBCFactory grins_factory_yz_symmetry("yz_symmetry");
  YZSymmetryBCFactory grins_factory_yz_symmetry_old_style("yz_symmetry_old_style");
  XZSymmetryBCFactory grins_factory_xz_symmetry("xz_symmetry");
  XZSymmetryBCFactory grins_factory_xz_symmetry_old_style("xz_symmetry_old_style");
  XYSymmetryBCFactory grins_factory_xy_symmetry("xy_symmetry");
  XYSymmetryBCFactory grins_factory_xy_symmetry_old_style("xy_symmetry_old_style");
  RollerXBCFactory grins_factory_roller_x("roller_x");
  RollerXBCFactory grins_factory_roller_x_old_style("roller_x_old_style");
  RollerYBCFactory grins_factory_roller_y("roller_y");
  RollerYBCFactory grins_factory_roller_y_old_style("roller_y_old_style");
  RollerZBCFactory grins_factory_roller_z("roller_z");
  RollerZBCFactory grins_factory_roller_z_old_style("roller_z_old_style");
  AxisymmetryBCFactory grins_factory_velocity_axisymmetry("Velocity_axisymmetric");
  AxisymmetryBCFactory grins_factory_velocity_axisymmetry_old_style("Velocity_axisymmetric_old_style");
} // end namespace GRINS
