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


// GRINS
#include "grins/collimated_gaussian_laser_intensity_profile.h"

namespace GRINS
{
  CollimatedGaussianLaserIntensityProfile::CollimatedGaussianLaserIntensityProfile(libMesh::Real w)
    : _w(w)
  {}

  void CollimatedGaussianLaserIntensityProfile::init( const std::vector<libMesh::Point> & quadrature_xyz,
                                                      const libMesh::Point & laser_centerline)
  {
    _intensity_vals.resize(quadrature_xyz.size());

    for (unsigned int i = 0; i < quadrature_xyz.size(); ++i)
      {
        libMesh::Real radius = ( quadrature_xyz[i] - laser_centerline).norm();
        _intensity_vals[i] = std::exp( -2.0 * (radius*radius)/(_w*_w) );;
      }

  }

}

