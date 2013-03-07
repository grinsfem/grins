//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id: bc_factory.C 33233 2012-09-21 05:22:22Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ERROR_ESTIMATOR_FACTORY_H
#define ERROR_ESTIMATOR_FACTORY_H

// GRINS
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/adjoint_refinement_estimator.h"

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
								 std::tr1::shared_ptr<QoIBase> qoi_base );

    virtual std::tr1::shared_ptr<libMesh::AdjointRefinementEstimator> build_adjref( const GetPot& input,
										    std::tr1::shared_ptr<QoIBase> qoi_base );

    private:

    libMesh::Real _refine_fraction;
    libMesh::Real _coarsen_fraction;

  };

} // end namespace GRINS
#endif // ERROR_ESTIMATOR_FACTORY_H
