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


#ifndef GRINS_INC_NAVIER_STOKES_BASE_H
#define GRINS_INC_NAVIER_STOKES_BASE_H

//GRINS
#include "grins/physics.h"
#include "grins/primitive_flow_fe_variables.h"


namespace GRINS
{

  //! Physics class for Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Incompressible Navier-Stokes equations.
   */
  class IncompressibleNavierStokesBase : public Physics
  {
  public:

    IncompressibleNavierStokesBase(const std::string& physics_name, const GetPot& input);

    ~IncompressibleNavierStokesBase();

    //! Initialization of Navier-Stokes variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

  protected:

    //! Physical dimension of problem
    /*! \todo Do we really need to cache this? */
    unsigned int _dim;

    PrimitiveFlowFEVariables _flow_vars;

    //! Material parameters, read from input
    /** \todo Create objects to allow for function specification */
    libMesh::Number _rho, _mu;
    
  private:
    IncompressibleNavierStokesBase();

  };

} //End namespace block

#endif // GRINS_INC_NAVIER_STOKES_BASE_H
