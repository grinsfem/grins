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


#ifndef GRINS_LASER_ABSORPTION_H
#define GRINS_LASER_ABSORPTION_H

// GRINS
#include "grins/multi_qoi_base.h"
#include "grins/fem_function_and_derivative_base.h"
#include "grins/spectroscopic_transmission.h"
#include "grins/laser_intensity_profile_base.h"

namespace GRINS
{
  /*!
    The LaserAbsorption class is a 2D/3D analogue of SpectroscopicAbsorption that
    can calculate the total absorbed laser intensity of a two-dimensional
    laser beam across a given flow field.
    
    It uses Gauss quadrature to integrate the initial and final laser intensity profiles.
    
    Each quadrature point has a corresponding SpectroscopicTransmission object to calculate the
    transmitted laser intensity along its respective optical path (i.e. RayfireMesh).
  */
  class LaserAbsorption : public MultiQoIBase
  {
  public:
    /*!
      Two-dimensional class constructor

      @param absorb An AbsorptionCoeff object that will be shared by all internal SpectroscopicTransmission classes
      @param top_origin A point on the mesh boundary at the top of the 2D laser
      @param centerline_origin A point on the mesh boundary at the centerline of the 2D laser
      @param bottom_origin A point on the mesh boundary at the bottom of the 2D laser
      @param theta angle in the xy-plane for the laser path, measured from the positive x-axis
                    with counterclockwise being positive; \f$-2\pi \leq \theta \leq +2\pi\f$
      @param n_quadrature_point The number of quadrature points used to integrate the intensity profile.
                                Each quadrature point will have a separate SpectroscopicTransmission class
      @param intensity_profile A LaserIntensityProfileBase object
      @param qoi_name The name of the QoI
    */
    LaserAbsorption(const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> & absorb,
                    const libMesh::Point & top_origin, const libMesh::Point & centerline_origin,
                    const libMesh::Point & bottom_origin,
                    libMesh::Real theta, unsigned int n_quadrature_points,
                    std::shared_ptr<LaserIntensityProfileBase> intensity_profile,   
                    const std::string & qoi_name);

      /*!
      Three-dimensional class constructor
      
      <b>top_origin, centerline_origin, and bottom_origin must not all be colinear</b>

      @param absorb An AbsorptionCoeff object that will be shared by all internal SpectroscopicTransmission classes
      @param top_origin A point on the mesh boundary on the outside edge of the circular 3D laser
      @param centerline_origin A point on the mesh boundary at the centerline of the 3D circular laser
      @param bottom_origin A different point on the mesh boundary on the outside edge of the circular 3D laser.
                          Must not be colinear with other 2 origin points
      @param theta angle in the xy-plane for the laser path, measured from the positive x-axis
                    with counterclockwise being positive; \f$-2\pi \leq \theta \leq +2\pi\f$
      @param phi angle from the positive z-axis, with \f$\phi = \pi/2\f$ being the xy-plane;
                  \f$ 0 \leq \phi \leq \pi \f$
      @param n_quadrature_point The number of quadrature points *in each dimension* used to integrate the intensity profile.
                                For example, n_quadrature_points='4' will result in a total of 16 quadrature points (4 in each dimension).
                                Each quadrature point will have a separate SpectroscopicTransmission class
      @param intensity_profile A LaserIntensityProfileBase object
      @param qoi_name The name of the QoI
    */
    LaserAbsorption(const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> & absorb,
                    const libMesh::Point & top_origin, const libMesh::Point & centerline_origin,
                    const libMesh::Point & bottom_origin,
                    libMesh::Real theta, libMesh::Real phi, unsigned int n_quadrature_points,
                    std::shared_ptr<LaserIntensityProfileBase> intensity_profile,   
                    const std::string & qoi_name);

    //! Just call the default copy constructor
    virtual QoIBase * clone() const;

    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index);

    //! AMR not yet supported
    virtual void element_qoi_derivative(AssemblyContext & context,
                                        const unsigned int qoi_index);

    //! Sum _qoi_vals from all processors and then do the gauss quadrature
    //! integration to find the total absorbed laser intensity
    virtual void parallel_op( const libMesh::Parallel::Communicator & communicator,
                              libMesh::Number & sys_qoi,
                              libMesh::Number & local_qoi );

  private:
    //! intensity profile object used to get the laser intensity at the quadrature points
    std::shared_ptr<LaserIntensityProfileBase> _intensity_profile;

    //! gauss quadrature weights cached here for use in parallel_op()
    std::vector<libMesh::Real> _quadrature_weights;

  };

}
#endif // GRINS_LASER_ABSORPTION_H

