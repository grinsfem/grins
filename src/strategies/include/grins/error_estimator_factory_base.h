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

#ifndef GRINS_ERROR_ESTIMATOR_FACTORY_BASE_H
#define GRINS_ERROR_ESTIMATOR_FACTORY_BASE_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/multiphysics_sys.h"
#include "grins/error_estimator_options.h"

// libMesh
#include "libmesh/error_estimator.h"

namespace GRINS
{
  // According to the standard, we need a declaration of the
  // specialization which precedes any automatic instantiation.
  template <> const GetPot* FactoryWithGetPot<libMesh::ErrorEstimator>::_input;

  //! Builds VariableBase objects
  /*! Most variable classes only require a GetPot object to the constructor, but
    others may require more information. Subclasses can dictate the necessary
    behavior. */
  class ErrorEstimatorFactoryBase : public FactoryWithGetPot<libMesh::ErrorEstimator>
  {
  public:
    ErrorEstimatorFactoryBase( const std::string& estimator_name )
      : FactoryWithGetPot<libMesh::ErrorEstimator>(estimator_name)
    {}

    ~ErrorEstimatorFactoryBase(){};

    static void set_system( MultiphysicsSystem& system )
    { _system = &system; }

    static void set_estimator_options( const ErrorEstimatorOptions& estimator_options )
    { _estimator_options = &estimator_options; }

  protected:

    //! Subclasses implement this method for building the ErrorEstimator object.
    virtual std::unique_ptr<libMesh::ErrorEstimator> build_error_estimator( const GetPot& input,
                                                                            MultiphysicsSystem& system,
                                                                            const ErrorEstimatorOptions& estimator_options ) =0;

    //! Cache pointer to system
    /*! We can't copy this so it must be a pointer. We do *not* own
      this so do not delete! */
    static MultiphysicsSystem* _system;

    //! Cache pointer to system
    /*! We do *not* own this so do not delete! */
    static const ErrorEstimatorOptions* _estimator_options;

  private:

    virtual std::unique_ptr<libMesh::ErrorEstimator> create();

  };

  inline
  std::unique_ptr<libMesh::ErrorEstimator> ErrorEstimatorFactoryBase::create()
  {
    if( !_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building ErrorEstimator!");
    if( !_system )
      libmesh_error_msg("ERROR: must call set_system() before building ErrorEstimator!");
    if( !_estimator_options )
      libmesh_error_msg("ERROR: must call set_estimator_options() before building ErrorEstimator!");

    std::unique_ptr<libMesh::ErrorEstimator>
      new_estimator = this->build_error_estimator( *_input, *_system, *_estimator_options );

    libmesh_assert(new_estimator);

    return new_estimator;
  }

} // end namespace GRINS

#endif // GRINS_ERROR_ESTIMATOR_FACTORY_BASE_H
