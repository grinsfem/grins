//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#ifndef AXISYM_MUSHY_ZONE_SOLIDIFICATION_H
#define AXISYM_MUSHY_ZONE_SOLIDIFICATION_H

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

  //! Adds axisymmetric mushy zone solidification source term
  /*!
    This class implements the axisymmetric mushy zone solidification to model a phase change
    from liquid to solid aimed at vacuum arc remelting (VAR) applications. The form is
    \f$ \mathbf{F} = -\mu ( \mathbf{u} - \mathbf{u}_{cast} )/K_{perm} \f$
    where
    \f$ \mu = \f$ viscosity,
    \f$ \mathbf{u} = \f$ is the velocity
    \f$ \mathbf{u}_{cast} = \f$ is the casting velocity
    \f$ K_{perm} = \f$ is the permeability of the medium.
    The permeability is an algebraic function of the liquid fraction, namely
    \f$ K_{perm} = (1-\phi_l)^2 A_{perm}/\phi_l^3 \f$
    where \f$ A_{perm} \f$ is a constant dependent on the characteristics of the medium. 
    The scalar phase function \f$ \phi_l = f(T_{melt}, \Delta T) \f$ where
    \f$ T_{melt} = \f$ the local melting temperature and 
    \f$ \Delta T = \f$ is the interface thickness.
   */
  class AxisymmetricMushyZoneSolidification : public Physics
  {
  public:
    
    AxisymmetricMushyZoneSolidification()
      : Physics()
    {};

    virtual ~AxisymmetricMushyZoneSolidification()
    {};

    //! Read options from GetPot input file.
    virtual void read_input_options( GetPot& input );

    //! Initialization of AxisymmetricMushyZoneSolidification variables
    /*!
      There are actually no extra variables
     */
    /*! \todo Perhaps put a default method in the base class so we don't have
      to overload with nothing? */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Register variables needed by AxisymmetricMushyZoneSolidification
    /*! This will register the temperature and velocity variables from
      the AxisymmetricIncompNavierStokes and AxisymmetricHeatTransfer classes.*/
    virtual void register_variable_indices(GRINS::VariableMap &global_map);

    // Context initialization
    /*! Doesn't do anything for AxisymmetricMushyZoneSolidification since there
      are no new variables created */
    virtual void init_context( libMesh::DiffContext &context );

    //! Source term contribution for AxisymmetricMushyZoneSolidification
    /*! This is the main part of the class. This will add the source term to
        the AxisymmetricIncompNavierStokes class.
     */
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricMushyZoneSolidification.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    //! No constraint terms for AxisymmetricMushyZoneSolidification.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricMushyZoneSolidification.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    //! No mass terms for AxisymmetricMushyZoneSolidification.
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system ); 

    //! No new variables, so no local map
    virtual void build_local_variable_map();

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    // Indices for each (registered/non-owned) variable;
    //! Index for registered r-velocity field
    RegtdVariableIndex _u_r_var;

    //! Index for registered z-velocity field
    RegtdVariableIndex _u_z_var;

    //! Index for registered temperature field
    RegtdVariableIndex _T_var;

    // Names of each registered variable in the system

    //! Name of registered r-velocity
    std::string _u_r_var_name;

    //! Name of registered z-velocity
    std::string _u_z_var_name;

    //! Name of registered temperature
    std::string _T_var_name;

    //! Function to compute liquid fraction \f$ \phi_l \f$
    /*! Computes the liquid fraction of the material. Generally,
        this function has the form \f$ \phi_l = f( T, T_{melt}, \Delta T ) \f$
	where \f$ T \f$ is the local temperature, \f$ T_{melt} \f$ is the given
	melting temperature of the material, and \f$ \Delta T \f$ is the
	desired "interface" thickness (given). In this implemenation,
	\f$ \phi_l \f$ is a Hermite cubic spline in the "melt zone". This
	choice was made to have continuous derivatives.
     */
    double compute_liquid_phi( const double T );
    
    //! Function to compute \f$ \frac{ d \phi_l}{dT} \f$
    double dphi_dT( const double T );

    //! Function to compute permeability \f$ K_{perm} \f$
    double compute_K_perm( const double T );

    //! Function to compute \f$ \frac{dK_{perm}}{dT} \f$
    double dKperm_dT( const double T );

    //! Melting temperature parameter \f$ T_{melt} \f$
    double _T_melt;

    //! "Interface" thickness \f$ \Delta T \f$
    /*! The liquid fraction takes the form \f$ \phi_l = f( T, T_{melt}, \Delta T ) \f$
        so that \f$ \phi_l = 0, T \le T_{melt} - \Delta T \f$ and 
	\f$ \phi_l = 1.0, T \ge T_{melt} + \Delta T \f$
     */
    double _delta_T;

    //! Permeability parameter
    double _A_perm;

    //! Numerical stability parameter for permeability
    double _eps;

    //! Viscosity
    /*! \todo Needs to be moved to a central location */
    double _mu;

    //! Casting velocity
    libMesh::Point _u_cast;

  }; // class AxisymmetricMushyZoneSolidification

} // namespace GRINS
#endif //AXISYM_MUSHY_ZONE_SOLIDIFICATION_H
