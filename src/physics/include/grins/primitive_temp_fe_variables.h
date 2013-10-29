//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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

#ifndef GRINS_PRIMITIVE_TEMP_FE_VARIABLES_H
#define GRINS_PRIMITIVE_TEMP_FE_VARIABLES_H

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

// GRINS
#include "grins/var_typedefs.h"
#include "grins/primitive_temp_variables.h"

namespace GRINS
{
  class PrimitiveTempFEVariables : public PrimitiveTempVariables
  {
  public:

    PrimitiveTempFEVariables( const GetPot& input, const std::string& physics_name );
    ~PrimitiveTempFEVariables();

    virtual void init( libMesh::FEMSystem* system );

  protected:

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family;

    //! Element orders, read from input
    libMeshEnums::Order _T_order;

  private:

    PrimitiveTempFEVariables();

  };

} // end namespace GRINS

#endif // GRINS_PRIMITIVE_TEMP_FE_VARIABLES_H
