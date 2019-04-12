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

#ifndef GRINS_ADJOINT_ERROR_ESTIMATOR_FACTORIES_H
#define GRINS_ADJOINT_ERROR_ESTIMATOR_FACTORIES_H

// GRINS
#include "grins/error_estimator_factory_base.h"

// libMesh
#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/adjoint_refinement_estimator.h"
#include "libmesh/qoi_set.h"

namespace GRINS
{
  template<typename EstimatorType>
  class AdjointErrorEstimatorFactoryBase : public ErrorEstimatorFactoryBase
  {
  public:
    AdjointErrorEstimatorFactoryBase( const std::string& estimator_name )
      : ErrorEstimatorFactoryBase(estimator_name)
    {}

    ~AdjointErrorEstimatorFactoryBase(){};

  protected:

    virtual std::unique_ptr<libMesh::ErrorEstimator>
    build_error_estimator( const GetPot& input, MultiphysicsSystem& system,
                           const ErrorEstimatorOptions& estimator_options )
    {
      std::unique_ptr<libMesh::ErrorEstimator>
        raw_error_estimator( new EstimatorType );

      // Now cast to derived type
      EstimatorType* adjoint_error_estimator =
        libMesh::cast_ptr<EstimatorType*>(raw_error_estimator.get() );

      // Set a QoISet
      libMesh::QoISet qoi_set(system);
      adjoint_error_estimator->qoi_set() = qoi_set;

      // Now set the options for the derived type
      this->set_adjoint_estimator_options( input, estimator_options, *adjoint_error_estimator );

      return raw_error_estimator;
    }

    virtual void set_adjoint_estimator_options( const GetPot& input,
                                                const ErrorEstimatorOptions& estimator_options,
                                                EstimatorType& estimator ) =0;
  };

  class AdjointResidualErrorEstimatorFactory : public AdjointErrorEstimatorFactoryBase<libMesh::AdjointResidualErrorEstimator>
  {
  public:

    AdjointResidualErrorEstimatorFactory( const std::string& estimator_name )
      : AdjointErrorEstimatorFactoryBase<libMesh::AdjointResidualErrorEstimator>(estimator_name)
    {}

    ~AdjointResidualErrorEstimatorFactory(){};

  protected:

    virtual void set_adjoint_estimator_options( const GetPot& input,
                                                const ErrorEstimatorOptions& estimator_options,
                                                libMesh::AdjointResidualErrorEstimator& estimator );
  };

  class AdjointRefinementErrorEstimatorFactory : public AdjointErrorEstimatorFactoryBase<libMesh::AdjointRefinementEstimator>
  {
  public:

    AdjointRefinementErrorEstimatorFactory( const std::string& estimator_name )
      : AdjointErrorEstimatorFactoryBase<libMesh::AdjointRefinementEstimator>(estimator_name)
    {}

    ~AdjointRefinementErrorEstimatorFactory(){};

  protected:

    virtual void set_adjoint_estimator_options( const GetPot& input,
                                                const ErrorEstimatorOptions& estimator_options,
                                                libMesh::AdjointRefinementEstimator& estimator );
  };
} // end namespace GRINS

#endif // GRINS_ADJOINT_ERROR_ESTIMATOR_FACTORIES_H
