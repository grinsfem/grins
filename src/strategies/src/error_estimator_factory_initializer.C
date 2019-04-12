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
#include "grins/error_estimator_factory_initializer.h"

// GRINS
#include "grins/error_estimator_factory_basic.h"
#include "grins/adjoint_error_estimator_factories.h"
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/kelly_error_estimator.h"

namespace GRINS
{
  ErrorEstimatorFactoryInitializer::ErrorEstimatorFactoryInitializer()
  {
    static ErrorEstimatorFactoryBasic<libMesh::KellyErrorEstimator>
      grins_factory_kelly_error_est(StrategiesParsing::kelly_error_estimator());

    static ErrorEstimatorFactoryBasic<libMesh::PatchRecoveryErrorEstimator>
      grins_factory_patch_recovery_error_est(StrategiesParsing::patch_recovery_error_estimator());

    static AdjointRefinementErrorEstimatorFactory
      grins_factory_adjoint_refinement_estimator(StrategiesParsing::adjoint_refinement_error_estimator());

    static AdjointResidualErrorEstimatorFactory
      grins_factory_adjoint_residual_estimator(StrategiesParsing::adjoint_residual_error_estimator());
  }
} // end namespace GRINS
