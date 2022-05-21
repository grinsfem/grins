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
#include "grins/composite_qoi.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/diff_context.h"

namespace GRINS
{
  CompositeQoI::CompositeQoI()
    : libMesh::DifferentiableQoI()
  {
    // We initialize these to false and then reset as needed by each QoI
    assemble_qoi_sides = false;
    assemble_qoi_elements = false;
  }

  CompositeQoI::CompositeQoI( const CompositeQoI & original )
    : libMesh::DifferentiableQoI(original)
  {
    _inactive_element_vars = original._inactive_element_vars;
    _inactive_side_vars = original._inactive_side_vars;

    for( auto & qoi : original._qois )
      {
        std::unique_ptr<QoIBase> clone(qoi->clone());
        this->add_qoi(std::move(clone));
      }
  }

  std::unique_ptr<libMesh::DifferentiableQoI> CompositeQoI::clone()
  {
    // We can't use std::make_unique here since we've made the
    // copy-constructor protected
    return std::unique_ptr<CompositeQoI>( new CompositeQoI(*this) );
  }

  void CompositeQoI::add_qoi( std::unique_ptr<QoIBase> qoi )
  {
    if( qoi->assemble_on_interior() )
      {
        this->assemble_qoi_elements = true;
      }

    if( qoi->assemble_on_sides() )
      {
        this->assemble_qoi_sides = true;
      }

    _qois.push_back(std::move(qoi));
  }

  void CompositeQoI::init_qoi_count( libMesh::System & sys )
  {
    sys.init_qois(_qois.size());
  }

  void CompositeQoI::init( const GetPot& input, const MultiphysicsSystem& system )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      _qois[q]->init(input,system,q);

    // Figure out which variables are active in our QoIs
    std::set<unsigned int> active_element_vars;
    std::set<unsigned int> active_side_vars;
    for( auto & qoi : _qois )
      qoi->register_active_vars(active_element_vars, active_side_vars);

    // Cache the inactive ones for reference during init_context
    std::vector<unsigned int> all_vars;
    system.get_all_variable_numbers(all_vars);

    if( this->assemble_qoi_elements )
      for( auto var : all_vars )
        if( active_element_vars.find(var) == active_element_vars.end() )
          _inactive_element_vars.insert(var);

    if( this->assemble_qoi_sides )
      for( auto var : all_vars )
        if( active_side_vars.find(var) == active_side_vars.end() )
          _inactive_side_vars.insert(var);
  }

  void CompositeQoI::init_context( libMesh::DiffContext& context )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( auto & qoi : _qois )
      qoi->init_context(c);

    if( this->assemble_qoi_elements )
      for( auto var : _inactive_element_vars )
        c.get_element_fe(var)->get_nothing();

    if( this->assemble_qoi_sides )
      for( auto var : _inactive_side_vars )
        c.get_side_fe(var)->get_nothing();
  }

  void CompositeQoI::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number>& param_pointer )
    const
  {
    for( auto & qoi : _qois )
      qoi->register_parameter(param_name, param_pointer);
  }

  void CompositeQoI::reinit(MultiphysicsSystem & system)
  {
    // call reinit() on each qoi
    for( auto & qoi : _qois )
      qoi->reinit(system);
  }

  void CompositeQoI::element_qoi( libMesh::DiffContext& context,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).element_qoi(c,q);
  }

  void CompositeQoI::element_qoi_derivative( libMesh::DiffContext& context,
                                             const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).element_qoi_derivative(c,q);
  }

  void CompositeQoI::side_qoi( libMesh::DiffContext& context,
                               const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).side_qoi(c,q);
  }

  void CompositeQoI::side_qoi_derivative( libMesh::DiffContext& context,
                                          const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).side_qoi_derivative(c,q);
  }

  void CompositeQoI::parallel_op( const libMesh::Parallel::Communicator& communicator,
                                  std::vector<libMesh::Number>& sys_qoi,
                                  std::vector<libMesh::Number>& local_qoi,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).parallel_op( communicator, sys_qoi[q], local_qoi[q] );
  }

  void CompositeQoI::thread_join( std::vector<libMesh::Number>& qoi,
                                  const std::vector<libMesh::Number>& other_qoi,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).thread_join( qoi[q], other_qoi[q] );
  }

  void CompositeQoI::finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives,
                                         std::size_t qoi_index)
  {
    for( auto & qoi : _qois )
      qoi->finalize_derivative(derivatives,qoi_index);
  }

  void CompositeQoI::output_qoi( std::ostream& out ) const
  {
    for( auto & qoi : _qois )
      qoi->output_qoi(out);
  }

  libMesh::Number CompositeQoI::get_qoi_value( unsigned int qoi_index ) const
  {
    return _qois[qoi_index]->value();
  }

} // end namespace GRINS
