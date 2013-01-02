//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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
#ifndef GRINS_CONSTANT_WITH_EXP_LAYER_H
#define GRINS_CONSTANT_WITH_EXP_LAYER_H

// libMesh
#include "fe_base.h"
#include "function_base.h"

namespace GRINS
{
  class ConstantWithExponentialLayer : public libMesh::FunctionBase<libMesh::Number>
  {
  public:

    ConstantWithExponentialLayer( const Real u0, const Real factor,
				  const Real r0, const Real delta );

    virtual ~ConstantWithExponentialLayer( );

    virtual libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > clone() const;

    virtual libMesh::Number operator()( const Point &p, const Real time );

    virtual void operator()( const Point &p, 
			     const Real time, 
			     libMesh::DenseVector<Number> &output );

    virtual libMesh::Number component( unsigned int i, const Point &p, 
				       const Real time );
    
  protected:
    
    Real eval( const Real r ) const;

    //! Coefficients defining parabola
    Real _u0;
    
    Real _factor;

    Real _r0;

    Real _delta;

  private:
    
    ConstantWithExponentialLayer();

  }; // class ConstantWithExponentialLayer
  
  inline
  libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> > ConstantWithExponentialLayer::clone() const
  {
    return libMesh::AutoPtr< libMesh::FunctionBase<libMesh::Number> >( new ConstantWithExponentialLayer( *this ) );
  }

  inline
  libMesh::Number ConstantWithExponentialLayer::operator()( const Point &p, 
							    const Real /*time*/ )
  {
    return this->eval( p(0) );
  }

  inline
  Real ConstantWithExponentialLayer::eval( const Real r ) const
  {
    return _u0*( 1.0 - std::exp( _factor*(r - _r0)/_delta) );
  }

  inline
  void ConstantWithExponentialLayer::operator()( const Point &p, 
						 const Real time, 
						 libMesh::DenseVector<libMesh::Number> &output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p, time);
      }
    return;
  }

  inline
  libMesh::Number ConstantWithExponentialLayer::component( unsigned int /*i*/,
							   const Point &p, 
							   const Real time )
  {
    return (*this)(p, time);
  }

  

} // namespace GRINS

#endif //GRINS_CONSTANT_WITH_EXP_LAYER_H
