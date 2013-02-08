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

// This class
#include "ignite_initial_guess.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"

namespace Bunsen
{

  template<class NumericType>
  IgniteInitialGuess<NumericType>::IgniteInitialGuess( const GetPot& input,
						       GRINS::MultiphysicsSystem& restart_system,
						       const GRINS::MultiphysicsSystem& init_system )
    : libMesh::FEMFunctionBase<NumericType>(),
      _restart_system( restart_system ),
      _T_var( init_system.variable_number( input( "VariableNames/temperature", "T") ) ),
      _r_min( input("InitialConditions/r_min", 0.0) ),
      _r_max( input("InitialConditions/r_max", 0.0) ),
      _z_min( input("InitialConditions/z_min", 0.0) ),
      _z_max( input("InitialConditions/z_max", 0.0) ),
      _T_value( input("InitialConditions/T_init", 0.0) )
  {
    /* ------ Build up variable map between systems ------ */
    std::vector<GRINS::VariableIndex> init_vars;

    // libMesh will resize the vector
    init_system.get_all_variable_numbers(init_vars);

    std::cout << "init_vars.size() = " << init_vars.size() << std::endl;
    
    for( unsigned int v = 0; v != init_vars.size(); v++ )
      {
	std::string var_name = init_system.variable_name( init_vars[v] );
	_var_map.insert( std::make_pair( init_vars[v], restart_system.variable_number(var_name) ) );
      }

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
    /*
    for( typename std::map<GRINS::VariableIndex,GRINS::VariableIndex>::const_iterator it = _var_map.begin();
	 it != _var_map.end(); it++ )
      {
	libMesh::FEBase* elem_fe = NULL;
	context.get_element_fe( it->second, elem_fe );
	elem_fe->get_phi();
	elem_fe->get_dphi();
	elem_fe->get_JxW();
	elem_fe->get_xyz();

	libMesh::FEBase* side_fe = NULL;
	context.get_side_fe( it->second, side_fe );
	side_fe->get_phi();
	side_fe->get_dphi();
	side_fe->get_JxW();
	side_fe->get_xyz();
      }
    */

    // Create the context we'll be using to compute MultiphysicsSystem quantities
    _restart_context.reset( new libMesh::FEMContext( _restart_system ) );
    _restart_system.init_context(*_restart_context);
    return;
  }

  template<class NumericType>
  NumericType IgniteInitialGuess<NumericType>::component( const libMesh::FEMContext& context, 
							  unsigned int component,
							  const libMesh::Point& p,
							  libMesh::Real /*time*/ )
  {
    // Check if the Elem is the same between the incoming context and the cached one.
    // If not, reinit the cached MultiphysicsSystem context
    if( &(context.get_elem()) != &(_restart_context->get_elem()) )
      {
	_restart_context->pre_fe_reinit(_restart_system,&context.get_elem());
	_restart_context->elem_fe_reinit();
      }

    const GRINS::VariableIndex restart_var = _var_map.find(component)->second;

    libMesh::Real value = 0.0;
    
    const libMesh::Real r = p(0);
    const libMesh::Real z = p(1);
    
    if( component == _T_var )
      {
	if( r >= _r_min &&
	    r <= _r_max &&
	    z >= _z_min && 
	    z <= _z_max )
	  {
	    value = _T_value;
	  }
	else
	  {
	    value = _restart_context->point_value( restart_var, p );
	  }
      }
    else
      {
	value = _restart_context->point_value( restart_var, p );
      }

    return value;
  }

  /* ------------------------- Instantiate -------------------------*/
  template class IgniteInitialGuess<libMesh::Real>;
} // namespace Bunsen
