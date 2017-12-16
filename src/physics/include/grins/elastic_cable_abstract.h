//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_CABLE_ABSTRACT_H
#define GRINS_ELASTIC_CABLE_ABSTRACT_H

//GRINS
#include "grins/solid_mechanics_abstract.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/fe_base.h"

namespace GRINS
{
  class ElasticCableAbstract : public SolidMechanicsAbstract
  {
  public:

    ElasticCableAbstract( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~ElasticCableAbstract(){};

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

  protected:

    //! Cross-sectional area of the cable
    libMesh::Real _A;

    //! Cable density
    libMesh::Real  _rho;

    const libMesh::FEGenericBase<libMesh::Real>* get_fe( const AssemblyContext& context );

  private:

    ElasticCableAbstract();

  };

  inline
  const libMesh::FEGenericBase<libMesh::Real>* ElasticCableAbstract::get_fe( const AssemblyContext& context )
  {
    // For this Physics, we need to make sure that we grab only the 1D elements
    return context.get_element_fe(_disp_vars.u(),1);
  }
}

#endif // GRINS_ELASTIC_CABLE_ABSTRACT_H
