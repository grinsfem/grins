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

#ifndef GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H
#define GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H

// GRINS
#include "grins_config.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/point.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMContext;
}

namespace GRINS
{  
  //! Adds Axisymmetric Boussinesq bouyancy source term
  /*!
    This class implements the Axisymmetric Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho_0 = \f$ reference density, 
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

    ~AxisymmetricBoussinesqBuoyancy();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization of AxisymmetricBoussinesqBuoyancy variables
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Source term contribution for AxisymmetricBoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
        the AxisymmetricIncompNavierStokes class.
     */
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _T_order, _V_order;

    // Indices for each variable;
    //! Index for r-velocity field
    VariableIndex _u_r_var;

    //! Index for z-velocity field
    VariableIndex _u_z_var;

    //! Index for temperature field
    VariableIndex _T_var;

    // Names of each variable in the system

    //! Name of r-velocity
    std::string _u_r_var_name;

    //! Name of z-velocity
    std::string _u_z_var_name;

    //! Name of temperature
    std::string _T_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _rho_ref;

    //! \f$ T_0 = \f$ reference temperature 
    libMesh::Number _T_ref;

    //! \f$ \beta_T = \f$ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    /* \todo This should be stashed in a singleton class and brought in from there */
    libMesh::Point _g;

  private:
    AxisymmetricBoussinesqBuoyancy();

  }; // class AxisymmetricBoussinesqBuoyancy

} // namespace GRINS
#endif // GRINS_AXISYM_BOUSSINESQ_BUOYANCY_H
