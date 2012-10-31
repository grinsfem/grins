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

#ifndef GRINS_VORTICITY_H
#define GRINS_VORTICITY_H

// GRINS
#include "qoi_base.h"
#include "variable_name_defaults.h"

namespace GRINS
{
  //! Vorticity QoI
  /*!
    This class implement a vorticity QoI that can be used to both compute
    QoI values and drive QoI-based adaptive refinement. Currently, this QoI
    is only implemented in 2D and will error if it detects a three-dimensional
    problem.
   */
  class Vorticity : public QoIBase
  {
  public:

    //! Constructor
    /*! Constructor takes GetPot object to read any input options associated
        with this QoI */
    Vorticity( const GetPot& input );

    virtual ~Vorticity();

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual libMesh::AutoPtr<libMesh::DifferentiableQoI> clone();

    //! Initialize local variables
    /*! Any local variables that need information from libMesh get initialized
        here. For example, variable indices. */
    virtual void init( const GetPot& input, const libMesh::FEMSystem& system );

    //! Compute the qoi value.
    /*! Currently, only implemented for 2D. Assumes that the vorticity will be
        computed over area of input subdomain id. Vorticity computed as 
        \f$ \int_{\Omega} \nabla \times \mathbf{u} \; d\mathbf{x}\f$*/
    virtual void element_qoi( DiffContext& context, const QoISet& qoi_indices );

    //! Compute the qoi derivative with respect to the solution.
    /*! Currently, only implemented for 2D. Assumes that the vorticity will be
        computed over area of input subdomain id. */
    virtual void element_qoi_derivative( DiffContext &context, const QoISet &qois );

  protected:

    virtual void read_input_options( const GetPot& input );

    //! u-velocity component variable index
    VariableIndex _u_var;

    //! v-velocity component variable index
    VariableIndex _v_var;

    //! List of sumdomain ids for which we want to compute this QoI
    std::set<libMesh::subdomain_id_type> _subdomain_ids;

  private:
    //! User never call default constructor.
    Vorticity();

  };
}
#endif //GRINS_VORTICITY_H
