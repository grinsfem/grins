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
// $Id:$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef AXISYM_ELECTROSTATICS_H
#define AXISYM_ELECTROSTATICS_H

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/parameters.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

//GRINS
#include "grins_config.h"
#include "grins/physics.h"
#include "axisym_electrostatics_bc_handling.h"

namespace GRINS
{

  //! Physics class for Axisymmetric Electrostatics
  /*
    This physics class implements electrostatics using the potential form of the
    equations.
  */
  class AxisymmetricElectrostatics : public Physics
  {
  public:

    AxisymmetricElectrostatics( const std::string& physics_name, const GetPot& input );
    ~AxisymmetricElectrostatics();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization  AxisymmetricHeatTransfer variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool request_jacobian,
					  libMesh::FEMContext& context );

    virtual void side_time_derivative( bool request_jacobian,
				       libMesh::FEMContext& context );

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    //! Index for electric potential
    VariableIndex _V_var;

    //! Name of r-velocity
    std::string _V_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _V_order;

  private:
    AxisymmetricElectrostatics();

  }; // class AxisymmetricElectrostatics
} //namespace GRINS

#endif //AXISYM_ELECTROSTATICS_H
