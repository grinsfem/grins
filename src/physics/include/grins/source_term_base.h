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


#ifndef GRINS_SOURCE_TERM_BASE_H
#define GRINS_SOURCE_TERM_BASE_H

// GRINS
#include "grins/physics.h"
#include "grins/grins_enums.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  //! Base class for generic source function term, f(x,y,z,t).
  /*! Idea is for the user to specify which variables they want to add a source
    term and then, depending on the subclass, parse for each variable and
    have it added to that equation.*/
  class SourceTermBase : public Physics
  {
  public:

    SourceTermBase( const std::string& physics_name, const GetPot& input );

    virtual ~SourceTermBase();

    virtual void init_variables( libMesh::FEMSystem* system );

  protected:

    std::vector<VariableIndex> _vars;
    std::vector<std::string> _var_names;
    std::vector<GRINSEnums::FEFamily> _var_FE;
    std::vector<GRINSEnums::Order> _var_order;

  private:

    SourceTermBase();

    //! Helper function
    void parse_var_info( const GetPot& input );

  };

} // end namespace GRINS

#endif // GRINS_SOURCE_TERM_BASE_H
