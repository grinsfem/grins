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

// This class
#include "grins/incompressible_hyperelasticity.h"

// GRINS
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"
#include "grins/cartesian_hyperelasticity.h"
#include "grins/incompressible_hyperelasticity_weak_form.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<unsigned int Dim,typename StrainEnergy>
  IncompressibleHyperelasticity<Dim,StrainEnergy>::IncompressibleHyperelasticity( const PhysicsName & physics_name,
                                                                                  const GetPot & input )
    :  HyperelasticityBase<Dim,StrainEnergy>(physics_name,physics_name,input),
    _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    _press_var.set_is_constraint_var(true);
    this->check_var_subdomain_consistency(_press_var);
  }

} // end namespace GRINS
