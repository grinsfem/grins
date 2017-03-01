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

// This class
#include "grins/solid_mechanics_abstract.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  SolidMechanicsAbstract::SolidMechanicsAbstract(const PhysicsName& physics_name,
                                                 const GetPot& input )
    : Physics(physics_name,input),
      _disp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<DisplacementVariable>(VariablesParsing::disp_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    // For solid mechanics problems, we need to set the sign for tractions
    // to '-' since the second order time solvers use a Newton residual of the form
    // M(u)\ddot{u} + C(u)\dot{u} + F(u) + G(u) = 0
    // In this case, then, the natural boundary conditions for the weak form
    // will all have a '-', so we need to set that so the input tractions are
    // positive.
    _disp_vars.set_neumann_bc_is_positive(false);

    this->check_var_subdomain_consistency(_disp_vars);
  }

  void SolidMechanicsAbstract::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march displacement forward in time, for as
    // many displacement variables as the dimension we're tracking

    system->time_evolving(_disp_vars.u(),2);

    if( this->_disp_vars.dim() > 1 )
      system->time_evolving(_disp_vars.v(),2);

    if( this->_disp_vars.dim() > 2 )
      system->time_evolving(_disp_vars.w(),2);
  }

} // end namespace GRINS
