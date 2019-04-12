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

// These classes
#include "grins/adjoint_error_estimator_factories.h"

// libMesh
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/enum_norm_type.h"

namespace GRINS
{
  void AdjointResidualErrorEstimatorFactory::set_adjoint_estimator_options( const GetPot& /*input*/,
                                                                            const ErrorEstimatorOptions& estimator_options,
                                                                            libMesh::AdjointResidualErrorEstimator& estimator )
  {
    libMesh::PatchRecoveryErrorEstimator*
      p1 = new libMesh::PatchRecoveryErrorEstimator;

    estimator.primal_error_estimator().reset( p1 );

    libMesh::PatchRecoveryErrorEstimator*
      p2 = new libMesh::PatchRecoveryErrorEstimator;

    estimator.dual_error_estimator().reset( p2 );

    bool patch_reuse = estimator_options.patch_reuse();
    estimator.primal_error_estimator()->error_norm.set_type( 0, libMesh::H1_SEMINORM );
    p1->set_patch_reuse( patch_reuse );

    estimator.dual_error_estimator()->error_norm.set_type( 0, libMesh::H1_SEMINORM );
    p2->set_patch_reuse( patch_reuse );
  }

  void AdjointRefinementErrorEstimatorFactory::set_adjoint_estimator_options( const GetPot& /*input*/,
                                                                              const ErrorEstimatorOptions& estimator_options,
                                                                              libMesh::AdjointRefinementEstimator& estimator )
  {
    estimator.number_h_refinements = estimator_options.n_adjoint_h_refinements();
    estimator.number_p_refinements = estimator_options.n_adjoint_p_refinements();
  }

} // end namespace GRINS
