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

#ifndef GRINS_BUNSEN_SOURCE_H
#define GRINS_BUNSEN_SOURCE_H

// GRINS
#include "grins/grins_enums.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

namespace GRINS
{
  class AssemblyContext;
}

namespace Bunsen
{
  class BunsenSource : public GRINS::Physics
  {
  public:

    BunsenSource( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~BunsenSource();

    virtual void read_input_options( const GetPot& input );

    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void element_time_derivative( bool compute_jacobian,
					  GRINS::AssemblyContext& context,
					  GRINS::CachedValues& cache );

  protected:

    const libMesh::Real _value;

    const libMesh::Real _r_max;
    const libMesh::Real _z_min;
    const libMesh::Real _z_max;

    GRINS::VariableIndex _T_var;
    std::string _T_var_name;

    GRINSEnums::FEFamily _T_FE_family;
    GRINSEnums::Order _T_order;

  private:

    BunsenSource();

  };

} // namespace Bunsen

#endif // GRINS_BUNSEN_SOURCE_H
