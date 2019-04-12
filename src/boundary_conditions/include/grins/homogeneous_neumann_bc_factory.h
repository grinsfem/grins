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

#ifndef GRINS_HOMOGENEOUS_NEUMANN_BC_FACTORY_H
#define GRINS_HOMOGENEOUS_NEUMANN_BC_FACTORY_H

// GRINS
#include "grins/neumann_bc_factory_abstract.h"

namespace GRINS
{
  //! Generic factory for homogeneous Neumann boundary conditions
  /*! Merely need to set _is_homogeneous to true at construction time, then
    the create() method in the base classes should just ignore this boundary
    condition factory and not construct anything, which is what we want
    for "do nothing" boundary conditions. */
  class HomogeneousNeumannBCFactory : public NeumannBCFactoryAbstract
  {
  public:
    HomogeneousNeumannBCFactory( const std::string& bc_type_name )
      : NeumannBCFactoryAbstract(bc_type_name)
    { _is_homogeneous = true; }

    virtual ~HomogeneousNeumannBCFactory(){};

  protected:

    //! Do nothing.
    /* This should never get called since we should be opt-ing out upstream of this
       function, so we error out if this gets called.*/
    virtual std::shared_ptr<NeumannBCAbstract>
    build_neumann_func( const GetPot& /*input*/,
                        MultiphysicsSystem& /*system*/,
                        const FEVariablesBase& /*fe_var*/,
                        const std::string& /*section*/ )
    { libmesh_error(); // This should never get called
      return std::shared_ptr<NeumannBCAbstract>(); }

  };

} // end namespace GRINS

#endif // GRINS_HOMOGENEOUS_NEUMANN_BC_FACTORY_H
