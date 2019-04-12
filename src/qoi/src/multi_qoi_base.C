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

// GRINS
#include "grins/multi_qoi_base.h"

namespace GRINS
{
  MultiQoIBase::MultiQoIBase(const std::string & qoi_name)
  : QoIBase(qoi_name)
  {
    _assemble_sides = false;
    _assemble_interior = false;
  }

  MultiQoIBase::~MultiQoIBase()
  {
    for(std::vector<std::unique_ptr<QoIBase>>::iterator qoi = _qois.begin(); qoi != _qois.end(); ++qoi)
        delete ((*qoi).release());

  }

  MultiQoIBase::MultiQoIBase(const MultiQoIBase & original)
    : QoIBase(original)
  {
    libmesh_assert_equal_to(original._qois.size(), original._qoi_vals.size());

    for (unsigned int q = 0; q < original._qois.size(); ++q)
      {
        this->add_qoi( *(original._qois[q].get()) );
        _qoi_vals[q] = original._qoi_vals[q];
      }

  }

  void MultiQoIBase::add_qoi(const QoIBase & qoi)
  {
    _qois.push_back(std::unique_ptr<QoIBase>(qoi.clone()));
    _qoi_vals.push_back(0.0);

    if(qoi.assemble_on_sides())
        _assemble_sides = true;

    if(qoi.assemble_on_interior())
        _assemble_interior = true;

  }

  QoIBase * MultiQoIBase::clone() const
  {
    MultiQoIBase * clone = new MultiQoIBase(_qoi_name);

    for(unsigned int q = 0; q < this->n_qois(); ++q)
        clone->add_qoi(this->get_qoi(q));

    return new MultiQoIBase(*clone);
  }

  void MultiQoIBase::init( const GetPot & input,
                     const MultiphysicsSystem & system,
                     unsigned int /*qoi_num*/)
  {
    for(unsigned int q = 0; q < this->n_qois(); ++q)
      _qois[q]->init(input,system,q);

  }

  void MultiQoIBase::init_context(AssemblyContext & context)
  {
    for(std::vector<std::unique_ptr<QoIBase>>::iterator qoi = _qois.begin(); qoi != _qois.end(); ++qoi)
        (*qoi)->init_context(context);

  }

  void MultiQoIBase::reinit(MultiphysicsSystem & system)
  {
    for (unsigned int q = 0; q < this->n_qois(); ++q)
      (this->get_qoi(q)).reinit(system);

  }

}

