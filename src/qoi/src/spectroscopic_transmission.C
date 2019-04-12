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
#include "grins/spectroscopic_transmission.h"
#include "grins/absorption_coeff.h"
#include "grins/integrated_function.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  SpectroscopicTransmission::SpectroscopicTransmission( const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real> > & absorb,
                                                    const std::shared_ptr<RayfireMesh> & rayfire, const std::string & qoi_name, bool output_as_csv)
    : SpectroscopicQoIBase(absorb,rayfire,qoi_name,output_as_csv)
  {}

  QoIBase * SpectroscopicTransmission::clone() const
  {
    return new SpectroscopicTransmission( *this );
  }

  void SpectroscopicTransmission::parallel_op( const libMesh::Parallel::Communicator & communicator,
                                             libMesh::Number & sys_qoi,
                                             libMesh::Number & local_qoi )
  {
    QoIBase::parallel_op(communicator,sys_qoi,local_qoi);

    // absorption coefficient is calculated in [cm^-1], but path length is given in [m]
    // 100.0 factor converts pathlength to [cm]
    sys_qoi = std::exp( -sys_qoi * 100.0 );
    QoIBase::_qoi_value = sys_qoi;
  }

  void SpectroscopicTransmission::finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives, std::size_t qoi_index)
  {
    if (!derivatives.closed())
      derivatives.close();

    // We recalculate the _qoi_value to make sure
    // it is set and up-to-date
    libMesh::QoISet qs;
    qs.add_index(qoi_index);
    _multiphysics_system->assemble_qoi(qs);

    derivatives.scale(100.0 * -QoIBase::_qoi_value);
  }

} //namespace GRINS
