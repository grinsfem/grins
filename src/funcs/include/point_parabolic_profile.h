//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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
#ifndef POINT_PARABOLIC_PROFILE_H
#define POINT_PARABOLIC_PROFILE_H

#include "point_func_base.h"

namespace GRINS
{
  
  //! Parabolic profile
  /*! Generic parabolic profile. Mainly used for defining inflow
      boundary conditions. Parabola takes the form:
  \f$ ax^2 + bxy + cy^2 + dx + ey + f \f$ */
  /** \todo Need to incorporate z-directions */
  class PointParabolicProfile : public BasePointFuncObj
  {
  public:
    
    //! Default constructor
    /*! Default constructor sets parameters for the profile:
      x_out = \f$ 4y(1-y) \f$
      y_out = 0
      z_out = 0 */
    PointParabolicProfile( );

    PointParabolicProfile( const double ax, const double bx, const double cx,
			   const double dx, const double ex, const double fx,
			   const double ay, const double by, const double cy,
			   const double dy, const double ey, const double fy,
			   const double az, const double bz, const double cz,
			   const double dz, const double ez, const double fz );

    virtual ~PointParabolicProfile( );

    virtual libMesh::Point operator()( const libMesh::Point& point );
    
  protected:
    
    inline double eval( const double a, const double b, const double c,
			const double d, const double e, const double f,
			const double x, const double y, const double z )
    {
      return a*x*x + b*x*y + c*y*y + d*x + e*y + f;
    };

    //! Coefficients defining parabola
    double _ax, _bx, _cx, _dx, _ex, _fx;
    double _ay, _by, _cy, _dy, _ey, _fy;
    double _az, _bz, _cz, _dz, _ez, _fz;

  }; // class PointParabolicProfile

} // namespace GRINS

#endif // POINT_PARABOLIC_PROFILE_H
