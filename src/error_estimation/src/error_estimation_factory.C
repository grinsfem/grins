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
// $Id: bc_factory.C 33233 2012-09-21 05:22:22Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/error_estimation_factory.h"

// libMesh
#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/getpot.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/qoi_set.h"

namespace GRINS
{

  ErrorEstimatorFactory::ErrorEstimatorFactory( )
  {
    return;
  }

  ErrorEstimatorFactory::~ErrorEstimatorFactory( )
  {
    return;
  }

  std::tr1::shared_ptr<ErrorEstimator> ErrorEstimatorFactory::build( const GetPot& input,
								     std::tr1::shared_ptr<GRINS::QoIBase> qoi_base )
  {
    // check if qoi_base is an empty pointer (user set no QoI), in that case return empty pointer.
    if( qoi_base == std::tr1::shared_ptr<QoIBase>() )
    {
      return std::tr1::shared_ptr<ErrorEstimator>();
    }

    std::tr1::shared_ptr<ErrorEstimator> error_estimator;
    AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;
    
    error_estimator.reset( adjoint_residual_estimator );
    
    //adjoint_residual_estimator->adjoint_already_solved = true;
    
    std::string estimator_type = input( "Adaptivity/estimator_type", "patch" );
    if( estimator_type == "patch" )
    {
      PatchRecoveryErrorEstimator *p1 = new PatchRecoveryErrorEstimator;
      adjoint_residual_estimator->primal_error_estimator().reset( p1 );
      
      PatchRecoveryErrorEstimator *p2 = new PatchRecoveryErrorEstimator;
      adjoint_residual_estimator->dual_error_estimator().reset( p2 );   
      
      bool patch_reuse = input( "Adaptivity/patch_reuse", true );
      adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type( 0, H1_SEMINORM );
      p1->set_patch_reuse( patch_reuse );
      
      adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type( 0, H1_SEMINORM );
      p2->set_patch_reuse( patch_reuse );
    }
    else
    {
      out << "Error: unrecognized option for estimator_type" << std::endl;
      libmesh_error();
    }
    
    return error_estimator;
  }

  std::tr1::shared_ptr<AdjointRefinementEstimator> ErrorEstimatorFactory::build_adjref(
    const GetPot& input,
    std::tr1::shared_ptr<GRINS::QoIBase> qoi_base )
  {
    std::tr1::shared_ptr<AdjointRefinementEstimator> error_estimator;

    /*
    error_estimator->_coarsen_fraction = input( "Adaptivity/coarsen_fraction", 0. );
    error_estimator->_refine_fraction = input( "Adaptivity/refine_fraction", 0. );

    AdjointRefinementEstimator *adjoint_refinement_estimator = new AdjointRefinementEstimator;
    
    error_estimator.reset( adjoint_refinement_estimator );
    
    adjoint_refinement_estimator->qoi_set() = qoi_base->get_enabled_qoi_set();
    // adjoint_refinement_estimator->adjoint_already_solved = true;

    // For now perform 1 uniform h-refinement as libMesh throws errors if adjoint space is same.
    adjoint_refinement_estimator->number_h_refinements = 1;
    */
        
    return error_estimator;
  }

} // namespace GRINS
