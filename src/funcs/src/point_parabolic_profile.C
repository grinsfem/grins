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

#include "point_parabolic_profile.h"

GRINS::PointParabolicProfile::PointParabolicProfile( )
  : BasePointFuncObj(),
    _ax(0.0), _bx(0.0), _cx(-4.0), _dx(0.0), _ex(4.0), _fx(0.0),
    _ay(0.0), _by(0.0), _cy(0.0), _dy(0.0), _ey(0.0), _fy(0.0),
    _az(0.0), _bz(0.0), _cz(0.0), _dz(0.0), _ez(0.0), _fz(0.0)
{
  return;
}

GRINS::PointParabolicProfile::PointParabolicProfile( const double ax, const double bx, const double cx,
						     const double dx, const double ex, const double fx,
						     const double ay, const double by, const double cy,
						     const double dy, const double ey, const double fy,
						     const double az, const double bz, const double cz,
						     const double dz, const double ez, const double fz )
  : BasePointFuncObj(),
    _ax(ax), _bx(bx), _cx(cx), _dx(dx), _ex(ex), _fx(fx),
    _ay(ay), _by(by), _cy(cy), _dy(dy), _ey(ey), _fy(fy),
    _az(az), _bz(bz), _cz(cz), _dz(dz), _ez(ez), _fz(fz)  
{
  return;
}

GRINS::PointParabolicProfile::~PointParabolicProfile()
{
  return;
}

libMesh::Point GRINS::PointParabolicProfile::operator()( const libMesh::Point& point )
{
  const double x = point(0);
  const double y = point(1);
  const double z = point(2);

  double x_out, y_out, z_out;

  x_out = this->eval( _ax, _bx, _cx, _dx, _ex, _fx, x, y, z );
  y_out = this->eval( _ay, _by, _cy, _dy, _ey, _fy, x, y, z );
  z_out = this->eval( _az, _bz, _cz, _dz, _ez, _fz, x, y, z );

  return libMesh::Point( x_out, y_out, z_out );
}
