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


#ifndef GRINS_LASER_INTENSITY_PROFILE_BASE_H
#define GRINS_LASER_INTENSITY_PROFILE_BASE_H

// C++
#include <vector>

// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/point.h"

namespace GRINS
{
  class LaserIntensityProfileBase
  {
  public:
    /*!
      This is a base class for handling evaluation of laser intensity profiles
      (used by the LaserAbsorption class)
    */

    virtual void init(const std::vector<libMesh::Point> & quadrature_xyz,
                      const libMesh::Point & laser_centerline) =0;

    libMesh::Real intensity(unsigned int index) const;

    virtual ~LaserIntensityProfileBase() = default;

  protected:
    std::vector<libMesh::Real> _intensity_vals;

  };

}
#endif // GRINS_LASER_INTENSITY_PROFILE_BASE_H
