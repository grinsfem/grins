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
#include "grins/spectroscopic_absorption.h"
#include "grins/absorption_coeff.h"
#include "grins/integrated_function.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/single_variable.h"
#include "grins/fem_function_and_derivative_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  SpectroscopicAbsorption::SpectroscopicAbsorption(const GetPot & input,const std::string & qoi_name,SharedPtr<FEMFunctionAndDerivativeBase<libMesh::Real> > absorb)
    : IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >(input,2 /* QGauss order */,absorb,"SpectroscopicAbsorption",qoi_name)
  {}

  QoIBase * SpectroscopicAbsorption::clone() const
  {
    return new SpectroscopicAbsorption( *this );
  }

  void SpectroscopicAbsorption::parallel_op( const libMesh::Parallel::Communicator & communicator,
                                             libMesh::Number & sys_qoi,
                                             libMesh::Number & local_qoi )
  {
    QoIBase::parallel_op(communicator,sys_qoi,local_qoi);

    // absorption coefficient is calculated in [cm^-1], but path length is given in [m]
    // 100.0 factor converts pathlength to [cm]
    sys_qoi = std::exp( -sys_qoi * 100.0 );
    QoIBase::_qoi_value = sys_qoi;
  }

} //namespace GRINS
