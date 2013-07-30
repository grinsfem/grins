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

#ifndef GRINS_ERROR_ESTIMATOR_FACTORY_H
#define GRINS_ERROR_ESTIMATOR_FACTORY_H

// GRINS
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/error_estimator.h"

// libMesh forward declartions
class GetPot;

namespace GRINS
{
  class ErrorEstimatorFactory
  {
  public:

    ErrorEstimatorFactory();

    virtual ~ErrorEstimatorFactory();

    virtual std::tr1::shared_ptr<libMesh::ErrorEstimator> build( const GetPot& input, 
                                                                 const libMesh::QoISet& qoi_set );

  protected:

    enum ErrorEstimatorEnum{ ADJOINT_RESIDUAL = 0,
                             ADJOINT_REFINEMENT,
                             KELLY,
                             PATCH_RECOVERY,
                             WEIGHTED_PATCH_RECOVERY,
                             UNIFORM_REFINEMENT };

    ErrorEstimatorEnum string_to_enum( const std::string& estimator_type ) const;

  };

} // end namespace GRINS
#endif // GRINS_ERROR_ESTIMATOR_FACTORY_H
