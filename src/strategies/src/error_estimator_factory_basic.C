//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/error_estimator_factory_basic.h"
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/kelly_error_estimator.h"

namespace GRINS
{
  template<typename EstimatorType>
  libMesh::UniquePtr<libMesh::ErrorEstimator>
  ErrorEstimatorFactoryBasic<EstimatorType>::build_error_estimator
  ( const GetPot& /*input*/, MultiphysicsSystem& /*system*/, const ErrorEstimatorOptions& /*estimator_options*/ )
  {
    return libMesh::UniquePtr<libMesh::ErrorEstimator>( new EstimatorType );
  }

  // Instantiate basic ErrorEstimator factories
  ErrorEstimatorFactoryBasic<libMesh::KellyErrorEstimator>
  grins_factory_kelly_error_est(StrategiesParsing::kelly_error_estimator());

  ErrorEstimatorFactoryBasic<libMesh::PatchRecoveryErrorEstimator>
  grins_factory_patch_recovery_error_est(StrategiesParsing::patch_recovery_error_estimator());

} // end namespace GRINS
