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


// This class
#include "grins/composite_qoi.h"

namespace GRINS
{
  CompositeQoI::CompositeQoI()
    : libMesh::DifferentiableQoI()
  {
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

  void CompositeQoI::add_qoi( QoIBase& qoi )
  {
    _qois.push_back( qoi.clone() );

    return;
  }

  void CompositeQoI::init_qoi( std::vector<Number>& sys_qoi )
  {
    sys_qoi.resize(qois.size(), 0.0);

    return;
  }

  void CompositeQoI::init( const GetPot& input, const MultiphysicsSystem& system )
  {
    for( std::vector<QoIBase*>::iterator qoi = _qois.begin();
         qoi != _qois.end(); ++qoi )
      {
        qoi->init(input,system);
      }

    return;
  }

  void CompositeQoI::parallel_op( const libMesh::Parallel::Communicator& communicator,
                                  std::vector<Number>& sys_qoi,
                                  std::vector<Number>& local_qoi,
                                  const QoISet& /*qoi_indices*/ )
  {
    for( unsigned int q = 0; q < _qois.size(); q++ )
      {
        (*_qois[q]).parallel_op( communicator, sys_qoi[q], local_qoi[q] );
      }

    return;
  }

  void CompositeQoI::output_qoi( std::ostream& out ) const
  {
    for( std::vector<QoIBase*>::iterator qoi = _qois.begin();
         qoi != _qois.end(); ++qoi )
      {
        qoi->output_qoi(out);
      }

    return;
  }

  libMesh::Number CompositeQoI::get_qoi( unsigned int qoi_index ) const
  {
    return (*_qois[qoi_index]).value();
  }

} // end namespace GRINS
