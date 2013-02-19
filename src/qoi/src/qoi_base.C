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
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  QoIBase::QoIBase()
    : libMesh::DifferentiableQoI()
  {
    return;
  }

  QoIBase::~QoIBase()
  {
    return;
  }

  void QoIBase::init_qoi( std::vector<Number>& sys_qoi )
  {
    sys_qoi.resize(1, 0.0);
    return;
  }

  void QoIBase::parallel_op( std::vector<Number>& sys_qoi, std::vector<Number>& local_qoi,
			     const QoISet& qoi_indices )
  {
    libMesh::DifferentiableQoI::parallel_op( sys_qoi, local_qoi, qoi_indices );
    _qoi_cache = sys_qoi;
    return;
  }

  void QoIBase::output_qoi( std::ostream& out ) const
  {
    if( !_qoi_cache.empty() )
      {
	out << "========================================================================" << std::endl;

	for(  unsigned int i = 0; i < _qoi_cache.size(); i++ )
	  {
	    out << "QoI #" << i << " = " 
		<< std::setprecision(16) 
		<< std::scientific
		<< _qoi_cache[i] << std::endl;
	  }

	out << "========================================================================" << std::endl;
      }

    return;
  }

  Number QoIBase::get_qoi( unsigned int qoi_index ) const
  {
    return _qoi_cache[qoi_index];
  }
  
} // namespace GRINS
