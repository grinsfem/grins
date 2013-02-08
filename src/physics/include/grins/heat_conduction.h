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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef GRINS_HEAT_CONDUCTION_H
#define GRINS_HEAT_CONDUCTION_H

#include "grins/physics.h"

namespace GRINS
{

  class HeatConduction : public Physics
  {

  public:

    HeatConduction( const GRINS::PhysicsName& physics_name, const GetPot& input );
    ~HeatConduction();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::FEMContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context,
				CachedValues& cache );

  protected:

    libMesh::Real forcing( const libMesh::Point& p );

    unsigned int _dim;

    //! Indices for each variable;
    VariableIndex _T_var;

    std::string _T_var_name;

    libMeshEnums::Order _T_order;

    libMeshEnums::FEFamily _T_FE_family;

    libMesh::Number _rho, _Cp, _k;

  private:

    HeatConduction();


  };

} // namespace GRINS

#endif // GRINS_HEAT_CONDUCTION_H
