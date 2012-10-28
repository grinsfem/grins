//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef AXISYM_LORENTZ_FORCE_H
#define AXISYM_LORENTZ_FORCE_H

#include "config.h"

#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"
#include "string_to_enum.h"
#include "fem_system.h"
#include "fem_context.h"

#include "physics.h"

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
  class AxisymmetricLorentzForce : public Physics
  {
  public:
    
    AxisymmetricLorentzForce( const std::string& physics_name, const GetPot& input );

    ~AxisymmetricLorentzForce();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization of AxisymmetricLorentzForce variables
    virtual void init_variables( libMesh::FEMSystem* system );

    // Context initialization
    /*! Doesn't do anything for AxisymmetricLorentzForce since there
      are no new variables registered */
    virtual void init_context( libMesh::DiffContext &context );

    //! Source term contribution for AxisymmetricLorentzForce
    /*! This is the main part of the class. This will add the source term to
        the AxisymmetricIncompNavierStokes class.
     */
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricLorentzForce.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    //! No constraint terms for AxisymmetricLorentzForce.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricLorentzForce.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    //! No mass terms for AxisymmetricLorentzForce.
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system ); 

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    //! Element type, read from input
    libMeshEnums::FEFamily _u_FE_family, _A_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _u_order, _A_order, _V_order;

    // Indices for each variable;
    //! Index for r-velocity field
    VariableIndex _u_r_var;

    //! Index for z-velocity field
    VariableIndex _u_z_var;

    //! Index for magnetic potential
    VariableIndex _A_var;

    //! Index for electric potential
    VariableIndex _V_var;

    // Names of each variable in the system

    //! Name of r-velocity
    std::string _u_r_var_name;

    //! Name of z-velocity
    std::string _u_z_var_name;

    //! Name of magnetic potential
    std::string _A_var_name;

    //! Name of electric potential
    std::string _V_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _sigma;

  private:
    AxisymmetricLorentzForce();

  }; // class AxisymmetricLorentzForce

} // namespace GRINS
#endif //AXISYM_LORENTZ_FORCE_H
