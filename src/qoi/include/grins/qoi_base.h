//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_QOI_BASE_H
#define GRINS_QOI_BASE_H

// C++
#include <iomanip>

// libMesh
#include "libmesh/diff_qoi.h"

// GRINS
#include "grins/parameter_user.h"
#include "grins/var_typedefs.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;
  class AssemblyContext;

  class QoIBase : public ParameterUser
  {
  public:

    QoIBase( const std::string& qoi_name );

    virtual ~QoIBase(){}

    //! Clone this QoI
    /*!
     * We return a raw pointer, but it is expected for the user to take ownership
     * and delete the object when done.
     */
    virtual QoIBase* clone() const =0;

    //! Does the QoI need an element interior assembly loop?
    /*! This is pure virtual to force to user to specify. */
    virtual bool assemble_on_interior() const =0;

    //! Does the QoI need a domain boundary assembly loop?
    /*! This is pure virtual to force to user to specify. */
    virtual bool assemble_on_sides() const =0;

    /*!
     * Method to allow QoI to cache any system information needed for QoI calculation,
     * for example, solution variable indices.
     */
    virtual void init( const GetPot& /*input*/,
                       const MultiphysicsSystem & /*system*/,
                       unsigned int /*qoi_num*/ ){}

    virtual void init_context( AssemblyContext& /*context*/ ){}

    //! Reinitialize QoI
    virtual void reinit(MultiphysicsSystem & /*system*/) {}

    //! Compute the qoi value for element interiors.
    /*! Override this method if your QoI is defined on element interiors */
    virtual void element_qoi( AssemblyContext& /*context*/,
                              const unsigned int /*qoi_index*/ ){}

    //! Compute the qoi derivative with respect to the solution on element interiors.
    /*! Override this method if your QoI is defined on element interiors */
    virtual void element_qoi_derivative( AssemblyContext & /*context*/,
                                         const unsigned int /*qoi_index*/ ){}

    //! Compute the qoi value on the domain boundary
    /*! Override this method if your QoI is defined on the domain boundary */
    virtual void side_qoi( AssemblyContext & /*context*/, const unsigned int /*qoi_index*/ ){}

    //! Compute the qoi derivative with respect to the solution on the domain boundary
    /*! Override this method if your QoI is defined on the domain boundary */
    virtual void side_qoi_derivative( AssemblyContext & /*context*/,
                                      const unsigned int /*qoi_index*/ ){}

    //! Call the parallel operation for this QoI and cache the value.
    /*!
     * By default, this is just a sum. Override if QoI is more complex.
     */
    virtual void parallel_op( const libMesh::Parallel::Communicator& communicator,
                              libMesh::Number& sys_qoi,
                              libMesh::Number& local_qoi );

    //! Call the operation to accumulate this QoI from multiple threads
    /*!
     * By default, this is just a sum. Override if QoI is more complex.
     */
    virtual void thread_join( libMesh::Number& qoi, const libMesh::Number& other_qoi );

    //! Finalize derivatives that require more than a simple sum.
    //! Does nothing by default.
    virtual void finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives);

    /*!
     * Basic output for computed QoI's. If fancier output is desired, override this method.
     */
    virtual void output_qoi( std::ostream& out ) const;

    //! Returns the current QoI value.
    libMesh::Number value() const;

    //! Returns the name of this QoI
    const std::string& name() const;

  protected:

    std::string _qoi_name;

    libMesh::Number _qoi_value;
  };

  inline
  libMesh::Number QoIBase::value() const
  {
    return _qoi_value;
  }

  inline
  const std::string& QoIBase::name() const
  {
    return _qoi_name;
  }

}
#endif // GRINS_QOI_BASE_H
