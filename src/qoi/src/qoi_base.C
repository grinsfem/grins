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
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

// GRINS
#include "grins/assembly_context.h"

namespace GRINS
{
  QoIBase::QoIBase( const std::string& qoi_name )
    : ParameterUser(qoi_name),
      _qoi_name(qoi_name),
      _qoi_value(0.0)
  {}

  void QoIBase::parallel_op( const libMesh::Parallel::Communicator& communicator,
                             libMesh::Number& sys_qoi,
                             libMesh::Number& local_qoi )
  {
    communicator.sum(local_qoi);

    sys_qoi = local_qoi;

    _qoi_value = sys_qoi;
  }

  void QoIBase::thread_join( libMesh::Number& qoi, const libMesh::Number& other_qoi )
  {
    qoi += other_qoi;
  }

  void QoIBase::finalize_derivative(libMesh::NumericVector<libMesh::Number> & /*derivatives*/)
  {
    // do nothing by default
  }

  void QoIBase::output_qoi( std::ostream& out ) const
  {
    out << _qoi_name+" = "
        << std::setprecision(16)
        << std::scientific
        << _qoi_value << std::endl;
  }

} // namespace GRINS
