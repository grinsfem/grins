//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

    return;
  }

  CompositeQoI::~CompositeQoI()
  {
    for( std::vector<QoIBase*>::iterator qoi = _qois.begin();
         qoi != _qois.end(); ++qoi )
      {
        delete (*qoi);
      }

    return;
  }

  libMesh::UniquePtr<libMesh::DifferentiableQoI> CompositeQoI::clone()
  {
    CompositeQoI* clone = new CompositeQoI;

    for( unsigned int q = 0; q < this->n_qois(); q++ )
      {
        clone->add_qoi( this->get_qoi(q) );
      }

    return libMesh::UniquePtr<libMesh::DifferentiableQoI>(clone);
  }

  void CompositeQoI::add_qoi( const QoIBase& qoi )
  {
    _qois.push_back( qoi.clone() );

    if( qoi.assemble_on_interior() )
      {
        this->assemble_qoi_elements = true;
      }

    if( qoi.assemble_on_sides() )
      {
        this->assemble_qoi_sides = true;
      }

    return;
  }

  void CompositeQoI::init_qoi( std::vector<libMesh::Number>& sys_qoi )
  {
    sys_qoi.resize(_qois.size(), 0.0);

    return;
  }

  void CompositeQoI::init( const GetPot& input, const MultiphysicsSystem& system )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      _qois[q]->init(input,system,q);
  }

  void CompositeQoI::init_context( libMesh::DiffContext& context )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( std::vector<QoIBase*>::iterator qoi = _qois.begin();
         qoi != _qois.end(); ++qoi )
      {
        (*qoi)->init_context(c);
      }

    return;
  }

  void CompositeQoI::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number>& param_pointer )
    const
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).register_parameter(param_name, param_pointer);
  }

  void CompositeQoI::reinit(MultiphysicsSystem & system)
  {
    // call reinit() on each qoi
    for (unsigned int i=0; i<this->n_qois(); i++)
      (this->get_qoi(i)).reinit(system);
  }

  void CompositeQoI::element_qoi( libMesh::DiffContext& context,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).element_qoi(c,q);
      }

    return;
  }

  void CompositeQoI::element_qoi_derivative( libMesh::DiffContext& context,
                                             const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).element_qoi_derivative(c,q);
      }

    return;
  }

  void CompositeQoI::side_qoi( libMesh::DiffContext& context,
                               const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).side_qoi(c,q);
      }

    return;
  }

  void CompositeQoI::side_qoi_derivative( libMesh::DiffContext& context,
                                          const libMesh::QoISet& /*qoi_indices*/ )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).side_qoi_derivative(c,q);
      }

    return;
  }

  void CompositeQoI::parallel_op( const libMesh::Parallel::Communicator& communicator,
                                  std::vector<libMesh::Number>& sys_qoi,
                                  std::vector<libMesh::Number>& local_qoi,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).parallel_op( communicator, sys_qoi[q], local_qoi[q] );
      }

    return;
  }

  void CompositeQoI::thread_join( std::vector<libMesh::Number>& qoi,
                                  const std::vector<libMesh::Number>& other_qoi,
                                  const libMesh::QoISet& /*qoi_indices*/ )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).thread_join( qoi[q], other_qoi[q] );
      }

    return;
  }

  void CompositeQoI::finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives, std::size_t /*qoi_index*/)
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      (*_qois[q]).finalize_derivative(derivatives);
  }

  void CompositeQoI::output_qoi( std::ostream& out ) const
  {
    for( std::vector<QoIBase*>::const_iterator qoi = _qois.begin();
         qoi != _qois.end(); ++qoi )
      {
        (*qoi)->output_qoi(out);
      }

    return;
  }

  libMesh::Number CompositeQoI::get_qoi_value( unsigned int qoi_index ) const
  {
    return (*_qois[qoi_index]).value();
  }

} // end namespace GRINS
