//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/postprocessed_quantities.h"

// GRINS
#include "grins/assembly_context.h"

namespace GRINS
{
  template<class NumericType>
  PostProcessedQuantities<NumericType>::PostProcessedQuantities( const GetPot& input )
    : libMesh::FEMFunctionBase<NumericType>()
  {
    if( input.have_variable("vis-options/output_vars") )
      {

        std::cerr << "================================================================================" << std::endl
                  << "WARNING: Detected input variable 'vis-options/output_vars." << std::endl
                  << "WARNING: THIS IS OUTDATED. output_vars is now input within" << std::endl
                  << "WARNING: each Physics section." << std::endl
                  << "         Erroring out to make you aware." << std::endl
                  << "================================================================================" << std::endl;
        libmesh_error();
      }
    return;
  }

  template<class NumericType>
  PostProcessedQuantities<NumericType>::~PostProcessedQuantities()
  {
    return;
  }

  template<class NumericType>
  unsigned int PostProcessedQuantities<NumericType>::register_quantity( std::string name )
  {
    // Check if this quantity has already been registered
    if( _quantity_name_index_map.find(name) != _quantity_name_index_map.end() )
      {
        std::cerr << "Error: trying to add existing quantity: " << name << std::endl;
        libmesh_error();
      }

    /* We add 1 so that the locally cached indices in each of the Physics
       can safely initialize the index to zero. In particular this ensures
       we count starting from 1. */
    unsigned int new_index = _quantity_name_index_map.size()+1;

    _quantity_name_index_map.insert( std::make_pair( name, new_index ) );

    return new_index;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::initialize( MultiphysicsSystem& system,
                                                         libMesh::EquationSystems& equation_systems )
  {
    // Only need to initialize if the user requested any output quantities.
    if( !_quantity_name_index_map.empty() )
      {
        // Need to cache the MultiphysicsSystem
        _multiphysics_sys = &system;

        libMesh::System& output_system = equation_systems.add_system<libMesh::System>("interior_output");

        for( std::map<std::string, unsigned int>::const_iterator it = _quantity_name_index_map.begin();
             it != _quantity_name_index_map.end();
             ++it )
          {
            unsigned int var = output_system.add_variable( it->first, libMesh::FIRST );
            _quantity_index_var_map.insert( std::make_pair(var,it->second) );
          }
      }

    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::update_quantities( libMesh::EquationSystems& equation_systems )
  {
    // Only do the projection if the user actually added any quantities to compute.
    if( !_quantity_name_index_map.empty() )
      {
        libMesh::System& output_system = equation_systems.get_system<libMesh::System>("interior_output");
        output_system.project_solution(this);
      }

    return;
  }


  template<class NumericType>
  NumericType PostProcessedQuantities<NumericType>::component( const libMesh::FEMContext& context,
                                                               unsigned int component,
                                                               const libMesh::Point& p,
                                                               libMesh::Real /*time*/ )
  {
    // Check if the Elem is the same between the incoming context and the cached one.
    // If not, reinit the cached MultiphysicsSystem context
    if(context.has_elem() && _multiphysics_context->has_elem())
      {
        if( &(context.get_elem()) != &(_multiphysics_context->get_elem()) )
          {
            _multiphysics_context->pre_fe_reinit(*_multiphysics_sys,&context.get_elem());
            _multiphysics_context->elem_fe_reinit();
          }
      }
    else if( !context.has_elem() && _multiphysics_context->has_elem() )
      {
        // Incoming context has NULL elem ==> SCALAR variables
        _multiphysics_context->pre_fe_reinit(*_multiphysics_sys,NULL);
        _multiphysics_context->elem_fe_reinit();
      }
    else if( context.has_elem() && !_multiphysics_context->has_elem() )
      {
        _multiphysics_context->pre_fe_reinit(*_multiphysics_sys,&context.get_elem());
        _multiphysics_context->elem_fe_reinit();
      }
    //else
    /* If has_elem() is false for both contexts, we're still dealing with SCALAR variables
       and therefore don't need to reinit. */

    libMesh::Real value = 0.0;

    // Quantity we want had better be there.
    libmesh_assert(_quantity_index_var_map.find(component) != _quantity_index_var_map.end());
    unsigned int quantity_index = _quantity_index_var_map.find(component)->second;

    _multiphysics_sys->compute_postprocessed_quantity( quantity_index,
                                                       *(this->_multiphysics_context),
                                                       p, value );

    return value;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::init_context( const libMesh::FEMContext& context )
  {
    // Make sure we prepare shape functions for our output variables.
    /*! \todo I believe this is redundant because it's done in the project_vector call. Double check. */
    for( typename std::map<VariableIndex,unsigned int>::const_iterator it = _quantity_index_var_map.begin();
         it != _quantity_index_var_map.end(); it++ )
      {
        libMesh::FEBase* elem_fe = NULL;
        context.get_element_fe( it->first, elem_fe );
        elem_fe->get_phi();
        elem_fe->get_dphi();
        elem_fe->get_JxW();
        elem_fe->get_xyz();

        libMesh::FEBase* side_fe = NULL;
        context.get_side_fe( it->first, side_fe );
        side_fe->get_phi();
        side_fe->get_dphi();
        side_fe->get_JxW();
        side_fe->get_xyz();
      }

    // Create the context we'll be using to compute MultiphysicsSystem quantities
    _multiphysics_context.reset( new AssemblyContext( *_multiphysics_sys ) );
    _multiphysics_sys->init_context(*_multiphysics_context);
    return;
  }

  template<class NumericType>
  void PostProcessedQuantities<NumericType>::operator()( const libMesh::FEMContext& context, const libMesh::Point& p,
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
  NumericType PostProcessedQuantities<NumericType>::operator()( const libMesh::FEMContext&,
                                                                const libMesh::Point&,
                                                                const libMesh::Real )
  {
    libmesh_error();
    return 0.0; //dummy
  }

  // Instantiate
  template class PostProcessedQuantities<libMesh::Real>;

} // namespace GRINS
