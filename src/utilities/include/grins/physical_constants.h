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


#ifndef GRINS_PHYSICAL_CONSTANTS_H
#define GRINS_PHYSICAL_CONSTANTS_H

#include "libmesh/libmesh_common.h"

namespace GRINS
{
  namespace Constants
  {
    /*!
     * Universal Gas Constant, R, expressed in J/(kmol-K)
     */
    const libMesh::Real R_universal = 8314.4621;

    /*!
     * Avogadro's number, particles per mole.
     */
    const libMesh::Real Avogadro = 6.02214179e23;

  } // namespace Constants
} // namespace GRINS

#endif //GRINS_PHYSICAL_CONSTANTS_H
