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

#ifndef GRINS_BOUSSINESQ_BUOYANCY_H
#define GRINS_BOUSSINESQ_BUOYANCY_H

// GRINS
#include "grins/physics.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/point.h"

namespace GRINS
{  
  //! Adds Boussinesq bouyancy source term
  /*!
    This class implements the Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho_0 = \f$ reference density, 
    \f$ T_0 = \f$ reference temperature,
    \f$ \beta_T = \f$ coefficient of thermal expansion, and
    \f$ \mathbf{g} = \f$ the gravitional vector.
    This source term to the governing flow equations through the
    element_time_derivative routine. This class requires a flow physics enabled
    and the ConvectiveHeatTransfer physics class enabled.
   */
  class BoussinesqBuoyancy : public Physics
  {
  public:
    
    BoussinesqBuoyancy( const std::string& physics_name, const GetPot& input );

    ~BoussinesqBuoyancy();

    //! Initialization of BoussinesqBuoyancy variables
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Source term contribution for BoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
        the IncompressibleNavierStokes class.
     */
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

  protected:

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _T_order, _V_order;

    //! Name of x-velocity
    std::string _u_var_name;

    //! Name of y-velocity
    std::string _v_var_name;

    //! Name of z-velocity
    std::string _w_var_name;

    //! Name of temperature
    std::string _T_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _rho_ref;

    //! \f$ T_0 = \f$ reference temperature 
    libMesh::Number _T_ref;

    //! \f$ \beta_T = \f$ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    libMesh::Point _g;

     //! Physical dimension of problem
    unsigned int _dim;

    // Indices for each variable;
    //! Index for x-velocity field
    VariableIndex _u_var;

    //! Index for y-velocity field
    VariableIndex _v_var;

    //! Index for z-velocity field
    VariableIndex _w_var;

    //! Index for temperature field
    VariableIndex _T_var;

  private:

    BoussinesqBuoyancy();

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_H
