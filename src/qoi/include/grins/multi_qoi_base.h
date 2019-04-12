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


#ifndef GRINS_MULTI_QOI_BASE_H
#define GRINS_MULTI_QOI_BASE_H

// GRINS
#include "grins/qoi_base.h"

namespace GRINS
{
  class MultiQoIBase : public QoIBase
  {
  public:

    /*!
      This class is an amalgamation of CompositeQoI and QoIBase.
      
      The intended purpose is to provide functionality for using multiple QoI objects
      to compute a single QoI value.
    */
    MultiQoIBase(const std::string & qoi_name);

    //! Release and delete all internally stored QoIBase objects
    virtual ~MultiQoIBase();

    //! Move all internal QoI object pointers by *deep copying* via clone() in add_qoi()
    MultiQoIBase(const MultiQoIBase & original);

    //! clone() and add QoI to internal vector
    void add_qoi(const QoIBase & qoi);

    unsigned int n_qois() const;

    const QoIBase & get_qoi(unsigned int qoi_index) const;

    QoIBase & get_qoi(unsigned int qoi_index);

    //! Uses clone() to deep copy all internal QoI objects
    virtual QoIBase * clone() const;

    virtual bool assemble_on_sides() const;

    virtual bool assemble_on_interior() const;

    //! init all internal QoI objects
    virtual void init( const GetPot & input,
                       const MultiphysicsSystem & system,
                       unsigned int qoi_num);

    //! call init_context on all internal QoI objects
    virtual void init_context(AssemblyContext & context);

    //! reinit all internal QoI objects
    virtual void reinit(MultiphysicsSystem & system);

  protected:
    //! A vector of internal QoIs that are *NOT* known to the context or CompositeQoI
    std::vector<std::unique_ptr<QoIBase>> _qois;

    //! Since the context does not know about the internal QoIs,
    //! we need to store their accumulated values separately
    std::vector<libMesh::Number> _qoi_vals;

    bool _assemble_sides;

    bool _assemble_interior;

  };

  inline
  bool MultiQoIBase::assemble_on_sides() const
  {
    return _assemble_sides;
  }

  inline
  bool MultiQoIBase::assemble_on_interior() const
  {
    return _assemble_interior;
  }

  inline
  unsigned int MultiQoIBase::n_qois() const
  {
    return _qois.size();
  }

  inline
  const QoIBase & MultiQoIBase::get_qoi(unsigned int qoi_index) const
  {
    libmesh_assert_less(qoi_index,this->n_qois());

    return (*(_qois[qoi_index].get()));
  }

  inline
  QoIBase & MultiQoIBase::get_qoi(unsigned int qoi_index)
  {
    libmesh_assert_less(qoi_index,this->n_qois());

    return (*(_qois[qoi_index]).get());
  }

}
#endif // GRINS_MULTI_QOI_BASE_H

