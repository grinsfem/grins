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
#include "grins/error_estimation_factory.h"

// libMesh
#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/getpot.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/qoi_set.h"
#include "libmesh/kelly_error_estimator.h"

namespace GRINS
{

  ErrorEstimatorFactory::ErrorEstimatorFactory()
  {
    return;
  }

  ErrorEstimatorFactory::~ErrorEstimatorFactory()
  {
    return;
  }

  std::tr1::shared_ptr<libMesh::ErrorEstimator> ErrorEstimatorFactory::build( const GetPot& input,
                                                                              const libMesh::QoISet& qoi_set )
  {
    std::string estimator_type = input("MeshAdaptivity/estimator_type", "patch_recovery");

    ErrorEstimatorEnum estimator_enum = this->string_to_enum( estimator_type );

    std::tr1::shared_ptr<libMesh::ErrorEstimator> error_estimator;

    switch( estimator_enum )
      {
      case(ADJOINT_RESIDUAL):
        {
          error_estimator.reset( new libMesh::AdjointResidualErrorEstimator );

          libMesh::AdjointResidualErrorEstimator*
            adjoint_error_estimator =
              libMesh::libmesh_cast_ptr<libMesh::AdjointResidualErrorEstimator*>
                ( error_estimator.get() );

          adjoint_error_estimator->qoi_set() = qoi_set;

          libMesh::PatchRecoveryErrorEstimator *p1 = new libMesh::PatchRecoveryErrorEstimator;
          adjoint_error_estimator->primal_error_estimator().reset( p1 );

          libMesh::PatchRecoveryErrorEstimator *p2 = new libMesh::PatchRecoveryErrorEstimator;
          adjoint_error_estimator->dual_error_estimator().reset( p2 );

          bool patch_reuse = input( "MeshAdaptivity/patch_reuse", false );
          adjoint_error_estimator->primal_error_estimator()->error_norm.set_type
            ( 0, libMesh::H1_SEMINORM );
          p1->set_patch_reuse( patch_reuse );

          adjoint_error_estimator->dual_error_estimator()->error_norm.set_type
            ( 0, libMesh::H1_SEMINORM );
          p2->set_patch_reuse( patch_reuse );
        }
        break;

      case(KELLY):
        {
          error_estimator.reset( new libMesh::KellyErrorEstimator );
        }
        break;

      case(PATCH_RECOVERY):
        {
          error_estimator.reset( new libMesh::PatchRecoveryErrorEstimator );
        }
        break;

      case(ADJOINT_REFINEMENT):
        {
          error_estimator.reset( new libMesh::AdjointRefinementEstimator );

          libMesh::AdjointRefinementEstimator*
            adjoint_error_estimator =
              libMesh::libmesh_cast_ptr<libMesh::AdjointRefinementEstimator*>
                ( error_estimator.get() );

          adjoint_error_estimator->qoi_set() = qoi_set;

          adjoint_error_estimator->number_h_refinements = input( "MeshAdaptivity/arefee_h_refs", 1 );
          adjoint_error_estimator->number_p_refinements = input( "MeshAdaptivity/arefee_h_refs", 0 );
        }
        break;

      case(WEIGHTED_PATCH_RECOVERY):
      case(UNIFORM_REFINEMENT):
        {
          libmesh_not_implemented();
        }
      break;

      // wat?!
      default:
        {
          libmesh_error();
        }

      } // switch( estimator_enum )

    return error_estimator;
  }

  ErrorEstimatorFactory::ErrorEstimatorEnum ErrorEstimatorFactory::string_to_enum( const std::string& estimator_type ) const
  {
    ErrorEstimatorEnum value;

    if( estimator_type == std::string("adjoint_residual") )
      {
        value = ADJOINT_RESIDUAL;
      }
    else if( estimator_type == std::string("adjoint_refinement") )
      {
        value = ADJOINT_REFINEMENT;
      }
    else if( estimator_type == std::string("kelly") )
      {
        value = KELLY;
      }
    else if( estimator_type == std::string("patch_recovery") )
      {
        value = PATCH_RECOVERY;
      }
    else if( estimator_type == std::string("weighted_patch_recovery") )
      {
        value = WEIGHTED_PATCH_RECOVERY;
      }
    else if( estimator_type == std::string("uniform_refinement") )
      {
        value = UNIFORM_REFINEMENT;
      }
    else
      {
        std::cerr << "Error: Invalid error estimator type " << estimator_type << std::endl
                  << "Valid error estimator types are: adjoint_residual" << std::endl
                  << "                                 adjoint_refinement" << std::endl
                  << "                                 kelly" << std::endl
                  << "                                 patch_recovery" << std::endl
                  << "                                 weighted_patch_recovery" << std::endl
                  << "                                 uniform_refinement" << std::endl;
        libmesh_error();
      }

    return value;
  }

} // namespace GRINS
