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


#ifndef GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H
#define GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H

// GRINS
#include "grins_config.h"
#include "grins/grins_enums.h"
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/point.h"

namespace GRINS
{
  //! Adds Axisymmetric Boussinesq bouyancy source term
  /*!
    This class implements the Axisymmetric Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho = \f$ density,
    \f$ T_0 = \f$ reference temperature,
    \f$ \beta_T = \f$ coefficient of thermal expansion, and
    \f$ \mathbf{g} = \f$ the gravitional vector.
    This source term is added to the governing flow equations through the
    element_time_derivative routine. This class requires an axisymmetric flow physics enabled
    and the AxisymmetricHeatTransfer physics class enabled.
  */
  class AxisymmetricBoussinesqBuoyancy : public Physics
  {
  public:

    AxisymmetricBoussinesqBuoyancy( const std::string& physics_name, const GetPot& input );

    ~AxisymmetricBoussinesqBuoyancy(){};

    virtual void init_context( AssemblyContext& context );

    //! Source term contribution for AxisymmetricBoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
      the AxisymmetricIncompNavierStokes class.
    */
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

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
    /* \todo This should be stashed in a singleton class and brought in from there */
    libMesh::Point _g;

  private:

    AxisymmetricBoussinesqBuoyancy();

    //! Read options from GetPot input file.
    void read_input_options( const GetPot& input );

  }; // class AxisymmetricBoussinesqBuoyancy

} // namespace GRINS
#endif // GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H
