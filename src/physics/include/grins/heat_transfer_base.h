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


#ifndef GRINS_HEAT_TRANSFER_BASE_H
#define GRINS_HEAT_TRANSFER_BASE_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"

namespace GRINS
{

  //! Physics class for Heat Transfer
  /*
    This physics class implements the classical Heat Transfer (neglecting viscous dissipation)
  */
  template<class Conductivity>
  class HeatTransferBase : public Physics
  {
  public:

    HeatTransferBase( const std::string& physics_name,
                      const std::string& core_physics_name,
                      const GetPot& input );

    ~HeatTransferBase(){};

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    VelocityVariable& _flow_vars;
    PressureFEVariable& _press_var;
    PrimitiveTempFEVariables& _temp_vars;

    //! Material parameters, read from input
    /*! \todo Need to generalize material parameters. Right now they
      are assumed constant */
    /*! \todo Shouldn't this rho be the same as the one in the flow? Need
      to figure out how to have those shared */
    libMesh::Number _rho, _Cp;

    //! Conductivity
    Conductivity _k;

  private:
    HeatTransferBase();
  };

} //End namespace block

#endif // GRINS_HEAT_TRANSFER_BASE_H
