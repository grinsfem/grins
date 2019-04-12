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

#ifndef GRINS_ERROR_ESTIMATOR_FACTORY_BASIC_H
#define GRINS_ERROR_ESTIMATOR_FACTORY_BASIC_H

#include "grins/error_estimator_factory_base.h"

namespace GRINS
{
  template<typename EstimatorType>
  class ErrorEstimatorFactoryBasic : public ErrorEstimatorFactoryBase
  {
  public:
    ErrorEstimatorFactoryBasic( const std::string& estimator_name )
      : ErrorEstimatorFactoryBase(estimator_name)
    {}

    ~ErrorEstimatorFactoryBasic(){};

  protected:

    virtual std::unique_ptr<libMesh::ErrorEstimator>
    build_error_estimator( const GetPot & /*input*/,
                           MultiphysicsSystem & /*system*/,
                           const ErrorEstimatorOptions & /*estimator_options*/ )
    {
      return std::unique_ptr<libMesh::ErrorEstimator>( new EstimatorType );
    }

  };

} // end namespace GRINS

#endif // GRINS_ERROR_ESTIMATOR_FACTORY_BASIC_H
