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

#ifndef GRINS_LORENTZ_FORCE_H
#define GRINS_LORENTZ_FORCE_H

// GRINS
#include "grins/physics.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
class GetPot;

namespace libMesh
{
  class FEMSystem;
  class FEMContext;
}

namespace GRINS
{
  class LorentzForce : public Physics
  {
  public:
    
    LorentzForce( const std::string& physics_name, const GetPot& input );

    ~LorentzForce();

    //! Initialization of LorentzForce variables
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Source term contribution for LorentzForce
    /*! This is the main part of the class. This will add the source term to
        the IncompNavierStokes class.
     */
    virtual void element_time_derivative( bool request_jacobian,
					  libMesh::FEMContext& context );

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    //! Element type, read from input
    libMeshEnums::FEFamily _u_FE_family, _A_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _u_order, _A_order, _V_order;

    // Indices for each variable;
    //! Index for x-velocity field
    VariableIndex _u_var;

    //! Index for y-velocity field
    VariableIndex _v_var;

    //! Index for z-velocity field
    VariableIndex _w_var;

    //! Index for magnetic potential
    VariableIndex _A_var;

    //! Index for electric potential
    VariableIndex _V_var;

    // Names of each variable in the system

    //! Name of x-velocity
    std::string _u_var_name;

    //! Name of y-velocity
    std::string _v_var_name;

    //! Name of z-velocity
    std::string _w_var_name;

    //! Name of magnetic potential
    std::string _A_var_name;

    //! Name of electric potential
    std::string _V_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _sigma;

  private:
    LorentzForce();

  }; // class LorentzForce

} // namespace GRINS

#endif // GRINS_LORENTZ_FORCE_H
