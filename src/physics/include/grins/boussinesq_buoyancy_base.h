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


#ifndef GRINS_BOUSSINESQ_BUOYANCY_BASE_H
#define GRINS_BOUSSINESQ_BUOYANCY_BASE_H

// GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/point.h"

namespace GRINS
{
  class BoussinesqBuoyancyBase : public Physics
  {
  public:

    BoussinesqBuoyancyBase( const std::string& physics_name, const GetPot& input );

    ~BoussinesqBuoyancyBase(){};

  protected:

    const VelocityVariable& _flow_vars;
    const PressureFEVariable& _press_var;
    const PrimitiveTempFEVariables& _temp_vars;

    //! \f$ \rho = \f$ density
    libMesh::Number _rho;

    //! \f$ T_0 = \f$ reference temperature
    libMesh::Number _T_ref;

    //! \f$ \beta_T = \f$ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    libMesh::Point _g;

  private:

    BoussinesqBuoyancyBase();

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_BASE_H
