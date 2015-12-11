//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_CABLE_BASE_H
#define GRINS_ELASTIC_CABLE_BASE_H

//GRINS
#include "grins/physics.h"
#include "grins/solid_mechanics_fe_variables.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/fe_base.h"

namespace GRINS
{
  class ElasticCableBase : public Physics
  {
  public:

    ElasticCableBase( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~ElasticCableBase();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

  protected:

    //! Cross-sectional area of the cable
    libMesh::Real _A;

    //! Cable density
    libMesh::Real  _rho;

    SolidMechanicsFEVariables _disp_vars;

    const libMesh::FEGenericBase<libMesh::Real>* get_fe( const AssemblyContext& context );

  private:

    ElasticCableBase();

  };

  inline
  const libMesh::FEGenericBase<libMesh::Real>* ElasticCableBase::get_fe( const AssemblyContext& context )
  {
    // For this Physics, we need to make sure that we grab only the 1D elements
    return context.get_element_fe(_disp_vars.u_var(),1);
  }
}

#endif // GRINS_ELASTIC_CABLE_BASE_H
