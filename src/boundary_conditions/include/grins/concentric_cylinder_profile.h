//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef CONCENTRIC_CYLINDER_PROFILE_H
#define CONCENTRIC_CYLINDER_PROFILE_H

// GRINS
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/function_base.h"

namespace GRINS
{
  
  //! Profile for flow between axially moving concentric cylinders
  /*! Profile generated between axially moving concentric cylinders.
      In this case, we assume the outer cylinder is stationary,
      and the inner cylinder is moving at speed u0. The profile
      is given by:
      \f$ u_0 \frac{ \log( r_1/r) }{\log( r_1/r_0)} \f$
      where: \f$ r_0 \f$ is the inner cylinder radius
      and \f$ r_1 \f$ is the outer cylinder radius. Note, that this
      assumes axisymmetry.
  */
  class ConcentricCylinderProfile : public libMesh::FunctionBase<libMesh::Number>
  {
  public:
    
    //! Default constructor
    /*! Default constructor sets parameters for the profile:
      value = \f$ 2.0 \frac{ \log( 2.0/r) }{\log(2.0)} \f$
    */
    ConcentricCylinderProfile( );

    ConcentricCylinderProfile( const double u0, const double r0, const double r1 );

    virtual ~ConcentricCylinderProfile( );

    virtual libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > clone() const;

    virtual libMesh::Number operator()( const libMesh::Point& p, const libMesh::Real time );

    virtual void operator()( const libMesh::Point& p, 
			     const libMesh::Real time, 
			     libMesh::DenseVector<libMesh::Number>& output );

    virtual libMesh::Number component( unsigned int i, const libMesh::Point& p, 
				       const libMesh::Real time );
    
  protected:
    
    inline 
    double eval( const double u0, const double r0, 
		 const double r1, const double r )
    {
      return u0*std::log(r1/r)/std::log(r1/r0);
    };

    //! Coefficients defining parabola
    double _u0, _r0, _r1;

  }; // class ConcentricCylinderProfile

} // namespace GRINS

#endif // CONCENTRIC_CYLINDER_H
