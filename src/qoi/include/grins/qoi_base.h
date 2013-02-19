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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_QOI_BASE_H
#define GRINS_QOI_BASE_H

// C++
#include <iomanip>

// libMesh
#include "libmesh/diff_qoi.h"

// GRINS
#include "grins/var_typedefs.h"

// libMesh forward declarations
class GetPot;

namespace libMesh
{
  class QoISet;
}

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  class QoIBase : public libMesh::DifferentiableQoI
  {
  public:

    QoIBase();
    virtual ~QoIBase();

    /*!
     * Method to allow QoI to cache any system information needed for QoI calculation,
     * for example, solution variable indices.
     */
    virtual void init( const GetPot& /*input*/, const MultiphysicsSystem& /*system*/ ){};

    /*!
     * Method to allow QoI to resize libMesh::System storage of QoI computations.
     * \todo Right now, we're only dealing with 1 QoI at a time. Need to generalize.
     */
    virtual void init_qoi( std::vector<libMesh::Number>& sys_qoi );

    /*!
     * We call the base class then grab the sys_qoi and cache it locally to output later.
     * If the QoI is not expressable as a sum over elements, then this will need to be
     * overridden with the correct libMesh::Parallel operations.
     */
    virtual void parallel_op( std::vector<libMesh::Number>& sys_qoi,
			      std::vector<libMesh::Number>& local_qoi,
			      const libMesh::QoISet& qoi_indices );

    /*!
     * Basic output for computed QoI's. If fancier output is desired, override this method.
     */
    virtual void output_qoi( std::ostream& out ) const;

    //! Accessor for value of QoI for given qoi_index.
    /*!
     * Returns value of QoI for qoi_index. Currently, we only store a single QoI,
     * so qoi_index should be zero.
     * \todo Maybe take a libMesh::QoISet instead?
     */
    libMesh::Number get_qoi( unsigned int qoi_index ) const;

  protected:

    std::vector<libMesh::Number> _qoi_cache;
  };
}
#endif // GRINS_QOI_BASE_H
