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
#ifndef BUNSEN_IGNITE_INITIAL_GUESS_H
#define BUNSEN_IGNITE_INITIAL_GUESS_H

// GRINS
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_context.h"

// GRINS forward declarations
namespace GRINS
{
  class MultiphysicsSystem;
}

// libmesh forward declarations
class GetPot;

namespace libMesh
{
  class Point;
}

namespace Bunsen
{

  template<class NumericType>
  class IgniteInitialGuess : public libMesh::FEMFunctionBase<NumericType>
  {
  public:
    IgniteInitialGuess( const GetPot& input,
		        GRINS::MultiphysicsSystem& restart_system,
			const GRINS::MultiphysicsSystem& init_system);

    virtual ~IgniteInitialGuess();

    /* Methods to override from FEMFunctionBase needed for libMesh-based evaluations */
    virtual void init_context( const libMesh::FEMContext & context);

    virtual libMesh::AutoPtr<libMesh::FEMFunctionBase<NumericType> > clone() const;

    virtual NumericType operator()( const libMesh::FEMContext& context, 
				    const libMesh::Point& p,
				    const libMesh::Real time = 0. );

    virtual void operator()( const libMesh::FEMContext& context, 
			     const libMesh::Point& p,
			     const libMesh::Real time,
			     libMesh::DenseVector<NumericType>& output );

    virtual NumericType component( const libMesh::FEMContext& context, 
				   unsigned int i,
				   const libMesh::Point& p,
				   libMesh::Real time=0. );

  protected:

    GRINS::MultiphysicsSystem& _restart_system;
    std::tr1::shared_ptr<libMesh::FEMContext> _restart_context;

    //! Map from init system variable number to restart system variable number
    std::map<GRINS::VariableIndex,GRINS::VariableIndex> _var_map;

    const GRINS::VariableIndex _T_var;

    const libMesh::Real _r_min;
    const libMesh::Real _r_max;
    const libMesh::Real _z_min;
    const libMesh::Real _z_max;

    const libMesh::Real _T_value;

  private:

    IgniteInitialGuess();

  };

  /* ------------------------- Inline Functions -------------------------*/

  template<class NumericType>
  inline
  libMesh::AutoPtr<libMesh::FEMFunctionBase<NumericType> > IgniteInitialGuess<NumericType>::clone() const
  {
    return libMesh::AutoPtr<libMesh::FEMFunctionBase<NumericType> >( new IgniteInitialGuess(*this) );
  }

  template<class NumericType>
  inline
  void IgniteInitialGuess<NumericType>::operator()( const libMesh::FEMContext& context,
						    const libMesh::Point& p,
						    const libMesh::Real time,
						    libMesh::DenseVector<NumericType>& output )
  {
    for( unsigned int i = 0; i != output.size(); i++ )
      {
	output(i) = this->component(context,i,p,time);
      }
    return;
  }

  template<class NumericType>
  inline
  NumericType IgniteInitialGuess<NumericType>::operator()( const libMesh::FEMContext&, 
							   const libMesh::Point&,
							   const libMesh::Real )
  {
    libmesh_error();
    return 0.0; //dummy
  }

} // namespace Bunsen

#endif // BUNSEN_IGNITE_INITIAL_GUESS_H
