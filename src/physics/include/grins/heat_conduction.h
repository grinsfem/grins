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

#ifndef GRINS_HEAT_CONDUCTION_H
#define GRINS_HEAT_CONDUCTION_H

// GRINS
#include "grins/physics.h"
#include "grins/single_variable.h"

namespace GRINS
{

  template<class Conductivity>
  class HeatConduction : public Physics
  {

  public:

    HeatConduction( const GRINS::PhysicsName& physics_name, const GetPot& input );
    ~HeatConduction(){};

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    PrimitiveTempFEVariables& _temp_vars;

    libMesh::Number _rho, _Cp;

    //! Conductivity
    Conductivity _k;

  private:

    HeatConduction();

  };

} // namespace GRINS

#endif // GRINS_HEAT_CONDUCTION_H
