//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef POINT_CONCENTRIC_CYLINDER_PROFILE_H
#define POINT_CONCENTRIC_CYLINDER_PROFILE_H

#include "point_func_base.h"

namespace GRINS
{
  
  //! Profile for flow between axially moving concentric cylinders
  /*! Profile generated between axially moving concentric cylinders.
      In this case, we assume the outer cylinder is stationary,
      and the inner cylinder is moving at speed u0. The profile
      is given by:
      \f$ u_0 \frac{ \ln( r_1/r) }{\ln( r_1/r_0} \f$
      where: \f$ r_0 \f$ is the inner cylinder radius
      and \f$ r_1 \f$ is the outer cylinder radius. Note, that this
      assumes axisymmetry.
  */
  class PointConcentricCylinderProfile : public BasePointFuncObj
  {
  public:
    
    //! Default constructor
    /*! Default constructor sets parameters for the profile:
      r_vel_out = 0.0
      z_out = \f$ 2.0 \frac{ \ln( 2.0/r) }{\ln(2.0)} \f$
    */
    PointConcentricCylinderProfile( );

    PointConcentricCylinderProfile( const double u0, 
				    const double r0, 
				    const double r1 );

    virtual ~PointConcentricCylinderProfile( );

    virtual libMesh::Point operator()( const libMesh::Point& point );
    
  protected:
    
    inline double eval( const double u0, const double r0, 
			const double r1, const double r )
    {
      return u0*std::log(r1/r)/std::log(r1/r0);
    };

    //! Coefficients defining parabola
    double _u0, _r0, _r1;

  }; // class PointConcentricCylinderProfile

} // namespace GRINS

#endif // POINT_CONCENTRIC_CYLINDER_H
