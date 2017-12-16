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


#ifndef GRINS_AXISYM_HEAT_TRANSFER_H
#define GRINS_AXISYM_HEAT_TRANSFER_H

//GRINS
#include "grins_config.h"
#include "grins/grins_enums.h"
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

namespace GRINS
{

  //! Physics class for Axisymmetric Heat Transfer
  /*
    This physics class implements the classical Axisymmetric Heat Transfer (neglecting viscous dissipation)
  */
  template<class Conductivity>
  class AxisymmetricHeatTransfer : public Physics
  {
  public:

    AxisymmetricHeatTransfer( const std::string& physics_name, const GetPot& input );

    ~AxisymmetricHeatTransfer(){};

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    // Registers all parameters in this physics and in its property
    // class
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    const VelocityVariable& _flow_vars;
    const PressureFEVariable& _press_var;
    const PrimitiveTempFEVariables& _temp_vars;

    //! Material parameters, read from input
    /*! \todo Need to generalize material parameters. Right now they
      are assumed constant */
    /*! \todo Shouldn't this rho be the same as the one in the flow? Need
      to figure out how to have those shared */
    libMesh::Number _rho, _Cp;

    Conductivity _k;

  private:

    AxisymmetricHeatTransfer();

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );

  };

} //End namespace block

#endif // GRINS_AXISYM_HEAT_TRANSFER_H
