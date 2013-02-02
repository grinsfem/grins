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
#ifndef PARABOLIC_PROFILE_H
#define PARABOLIC_PROFILE_H

// GRINS
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/function_base.h"

namespace GRINS
{
  
  //! Parabolic profile
  /*! Generic parabolic profile. Mainly used for defining inflow
      boundary conditions. Parabola takes the form:
  \f$ ax^2 + bxy + cy^2 + dx + ey + f \f$ */
  /** \todo Need to incorporate z-directions */
  class ParabolicProfile : public libMesh::FunctionBase<libMesh::Number>
  {
  public:
    
    //! Default constructor
    /*! Default constructor sets parameters for the profile:
      value = \f$ 4y(1-y) \f$ */
    ParabolicProfile( );

    ParabolicProfile( const double a, const double b, const double c,
		      const double d, const double e, const double f );

    virtual ~ParabolicProfile( );

    virtual libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > clone() const;

    virtual libMesh::Number operator()( const libMesh::Point &p, const libMesh::Real time );

    virtual void operator()( const libMesh::Point &p, 
			     const libMesh::Real time, 
			     libMesh::DenseVector<libMesh::Number> &output );

    virtual libMesh::Number component( unsigned int i, const libMesh::Point &p, 
				       const libMesh::Real time );
    
  protected:
    
    inline double eval( const double a, const double b, const double c,
			const double d, const double e, const double f,
			const double x, const double y )
    {
      return a*x*x + b*x*y + c*y*y + d*x + e*y + f;
    };

    //! Coefficients defining parabola
    double _a, _b, _c, _d, _e, _f;

  }; // class ParabolicProfile

} // namespace GRINS

#endif // PARABOLIC_PROFILE_H
