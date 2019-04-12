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


#ifndef GRINS_COLLIMATED_GAUSSIAN_LASER_INTENSITY_PROFILE_H
#define GRINS_COLLIMATED_GAUSSIAN_LASER_INTENSITY_PROFILE_H

// GRINS
#include "grins/laser_intensity_profile_base.h"

namespace GRINS
{
  class CollimatedGaussianLaserIntensityProfile : public LaserIntensityProfileBase
  {
  public:
    /*!
      Collimated gaussian laser intensity profile where the "spot size", w, is constant
      along the optical path. The 'spot size" is the radius at which the laser intensity
      drops to ~13.5% of its centerline value

      \f$ I(r) \ = \ I(0) \exp(-2 \frac{r^2}{w^2} ) \f$
    */
    CollimatedGaussianLaserIntensityProfile(libMesh::Real w);

    virtual void init(const std::vector<libMesh::Point> & quadrature_xyz,
                      const libMesh::Point & laser_centerline);

  private:
    libMesh::Real _w;
  };

}
#endif // GRINS_COLLIMATED_GAUSSIAN_LASER_INTENSITY_PROFILE_H

