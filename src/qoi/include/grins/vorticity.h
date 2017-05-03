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


#ifndef GRINS_VORTICITY_H
#define GRINS_VORTICITY_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  class VelocityVariable;

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
    Vorticity( const std::string& qoi_name );

    virtual ~Vorticity();

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    //! Initialize local variables
    /*! Any local variables that need information from libMesh get initialized
      here. For example, variable indices. */
    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

    virtual void init_context( AssemblyContext& context );

    //! Compute the qoi value.
    /*! Currently, only implemented for 2D. Assumes that the vorticity will be
      computed over area of input subdomain id. Vorticity computed as
      \f$ \int_{\Omega} \nabla \times \mathbf{u} \; d\mathbf{x}\f$*/
    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index );

    //! Compute the qoi derivative with respect to the solution.
    /*! Currently, only implemented for 2D. Assumes that the vorticity will be
      computed over area of input subdomain id. */
    virtual void element_qoi_derivative( AssemblyContext& context,
                                         const unsigned int qoi_index );

  protected:

    const VelocityVariable * _flow_vars;

    //! List of sumdomain ids for which we want to compute this QoI
    std::set<libMesh::subdomain_id_type> _subdomain_ids;

  private:
    //! User never call default constructor.
    Vorticity();

  };

  inline
  bool Vorticity::assemble_on_interior() const
  {
    return true;
  }

  inline
  bool Vorticity::assemble_on_sides() const
  {
    return false;
  }
}
#endif //GRINS_VORTICITY_H
