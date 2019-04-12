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

#ifndef GRINS_COMPOSITE_QOI_H
#define GRINS_COMPOSITE_QOI_H

// C++
#include <vector>
#include <ostream>

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/diff_qoi.h"
#include "libmesh/auto_ptr.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class DiffContext;
  class QoISet;
  namespace Parallel
  {
    class Communicator;
  }
}

// GRINS
#include "grins/qoi_base.h"

namespace GRINS
{
  // GRINS forward declarations
  class MultiphysicsSystem;

  class CompositeQoI : public libMesh::DifferentiableQoI
  {
  public:
    CompositeQoI();

    virtual ~CompositeQoI();

    //! Required to provide clone for adding QoI object to libMesh objects.
    /*! Note that we do a deep copy here since the previous object might
      get destroyed and wipe out the objects being pointed to in _qois. */
    virtual std::unique_ptr<libMesh::DifferentiableQoI> clone();

    virtual void add_qoi( const QoIBase& qoi );

    unsigned int n_qois() const;

    //! Each QoI will register its copy(s) of an independent variable
    //  named in this call.
    void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number>& param_pointer)
      const;

    /*!
     * Method to allow QoI to cache any system information needed for QoI calculation,
     * for example, solution variable indices.
     */
    virtual void init( const GetPot& input, const MultiphysicsSystem& system );

    /*!
     * Method to allow QoI to resize libMesh::System storage of QoI computations.
     */
    virtual void init_qoi( std::vector<libMesh::Number>& sys_qoi );

    virtual void init_context( libMesh::DiffContext& context );

    //! Reinitialize qoi
    virtual void reinit(MultiphysicsSystem & system);

    //! Compute the qoi value for element interiors.
    virtual void element_qoi( libMesh::DiffContext& context,
                              const libMesh::QoISet& qoi_indices );

    //! Compute the qoi derivative with respect to the solution on element interiors.
    virtual void element_qoi_derivative( libMesh::DiffContext &context,
                                         const libMesh::QoISet &qoi_indices );

    //! Compute the qoi value on the domain boundary
    virtual void side_qoi( libMesh::DiffContext& context, const libMesh::QoISet& qoi_indices );

    //! Compute the qoi derivative with respect to the solution on the domain boundary
    virtual void side_qoi_derivative( libMesh::DiffContext &context, const libMesh::QoISet &qois );

    //! Operation to accumulate the QoI from multiple MPI processes.
    /*!
     * Calls each QoI's parallel_op function.
     */
    virtual void parallel_op( const libMesh::Parallel::Communicator& communicator,
                              std::vector<libMesh::Number>& sys_qoi,
                              std::vector<libMesh::Number>& local_qoi,
                              const libMesh::QoISet& qoi_indices );

    //! Operation to accumulate the QoI from multiple MPI processes.
    /*!
     * Calls each QoI's thread_join function.
     */
    virtual void thread_join( std::vector<libMesh::Number>& qoi,
                              const std::vector<libMesh::Number>& other_qoi,
                              const libMesh::QoISet& qoi_indices );

    // Calls each QoI's finalize_derivative function
    virtual void finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives, std::size_t qoi_index);

    //! Basic output for computed QoI's.
    void output_qoi( std::ostream& out ) const;

    //! Accessor for value of QoI for given qoi_index.
    libMesh::Number get_qoi_value( unsigned int qoi_index ) const;

    const QoIBase& get_qoi( unsigned int qoi_index ) const;

    //! Non-const version needed for reinit()
    QoIBase& get_qoi( unsigned int qoi_index );

  protected:

    std::vector<QoIBase*> _qois;

  };

  inline
  unsigned int CompositeQoI::n_qois() const
  {
    return _qois.size();
  }

  inline
  const QoIBase& CompositeQoI::get_qoi( unsigned int qoi_index ) const
  {
    libmesh_assert_less( qoi_index, this->n_qois() );

    return (*_qois[qoi_index]);
  }

  inline
  QoIBase& CompositeQoI::get_qoi( unsigned int qoi_index )
  {
    libmesh_assert_less( qoi_index, this->n_qois() );

    return (*_qois[qoi_index]);
  }

} // end namespace GRINS

#endif // GRINS_COMPOSITE_QOI_H
