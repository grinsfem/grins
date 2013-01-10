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

#include "ignite_initial_guess.h"

namespace Bunsen
{

  template<class NumericType>
  IgniteInitialGuess<NumericType>::IgniteInitialGuess( const GetPot& input,
						       GRINS::MultiphysicsSystem& restart_system,
						       const GRINS::MultiphysicsSystem& init_system )
    : _restart_system( restart_system ),
      _prev_point(1.0e15,1.0e15,1.0e15) //Initialize to an absurd value
  {
    return;
  }

  template<class NumericType>
  IgniteInitialGuess<NumericType>::~IgniteInitialGuess()
  {
    return;
  }

  template<class NumericType>
  void IgniteInitialGuess<NumericType>::init_context( const libMesh::FEMContext& context )
  {
    // Create the context we'll be using to compute MultiphysicsSystem quantities
    _restart_context.reset( new libMesh::FEMContext( *_multiphysics_sys ) );
    _restart_sys->init_context(*_restart_context);
    return;
  }

  template<class NumericType>
  NumericType IgniteInitialGuess<NumericType>::component( const libMesh::FEMContext& context, 
							  unsigned int component,
							  const libMesh::Point& p,
							  Real /*time*/ )
  {
    // Check if the Elem is the same between the incoming context and the cached one.
    // If not, reinit the cached MultiphysicsSystem context
    if( &(context.get_elem()) != &(_restart_context->get_elem()) )
      {
	_restart_context->pre_fe_reinit(*_multiphysics_sys,&context.get_elem());
	_restart_context->elem_fe_reinit();
      }

    /* Optimization since we expect this function to be called many times with
       the same point. _prev_point initialized to something absurd so this should 
       always be false the first time. */
    if( _prev_point != p )
      {
	_prev_point = p;
	std::vector<libMesh::Point> point_vec(1,p);
	this->_cache.clear();
	_multiphysics_sys->compute_element_cache( *(this->_restart_context), point_vec, this->_cache );
      }

    const unsigned int restart_var = _var_map.find(component)->second;

    return _restart_context->point_value( restart_var, p );
  }

  /* ------------------------- Instantiate -------------------------*/
  template class IgniteInitialGuess<Real>;
} // namespace Bunsen
