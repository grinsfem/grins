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

#ifndef GRINS_THERMO_PRESSURE_FE_VARIABLE_H
#define GRINS_THERMO_PRESSURE_FE_VARIABLE_H

// GRINS
#include "grins/single_fe_type_variable.h"
#include "grins/thermo_pressure_variable.h"

namespace GRINS
{

  class ThermoPressureFEVariable : public SingleFETypeVariable,
                                   public ThermoPressureVariable
  {
  public:

    ThermoPressureFEVariable( const GetPot& input, const std::string& physics_name )
      :  SingleFETypeVariable(input,physics_name,"",this->subsection(),"SCALAR","FIRST"),
         ThermoPressureVariable(input)
    {
      // Currently only support SCALAR and FIRST
      libmesh_assert_equal_to( _family[0], libMesh::SCALAR );
      libmesh_assert_equal_to( _order[0], libMesh::FIRST );
    }

    ~ThermoPressureFEVariable(){};

    virtual void init( libMesh::FEMSystem* system )
    { this->default_fe_init(system, _var_names, _vars ); }

  private:

    ThermoPressureFEVariable();

  };

} // end namespace GRINS

#endif // GRINS_THERMO_PRESSURE_FE_VARIABLE_H
