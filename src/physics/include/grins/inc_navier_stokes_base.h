//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/velocity_fe_variables.h"
#include "grins/pressure_fe_variable.h"

namespace GRINS
{

  //! Physics class for Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Incompressible Navier-Stokes equations.
    This is a templated class, the class Viscosity can be instantiated as a specific type
    (right now:ConstantViscosity or SpatiallyVaryingViscosity) to allow the user
    to specify a constant or spatially varying viscosity in the input file
   */
  template<class Viscosity>
  class IncompressibleNavierStokesBase : public Physics
  {
  public:

    IncompressibleNavierStokesBase(const std::string& my_physics_name,
                                   const std::string& core_physics_name,
                                   const GetPot& input);

    ~IncompressibleNavierStokesBase(){};

    //! Initialization of Navier-Stokes variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

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

    // A getter function for the Viscosity object
    libMesh::Real get_viscosity_value(AssemblyContext& context, unsigned int qp) const;

  protected:

    //! Physical dimension of problem
    /*! \todo Do we really need to cache this? */
    unsigned int _dim;

    VelocityFEVariables _flow_vars;
    PressureFEVariable _press_var;
    //! Material parameters, read from input
    /** \todo Create objects to allow for function specification */
    libMesh::Number _rho;

    //! Viscosity object
    Viscosity _mu;

  private:
    IncompressibleNavierStokesBase();

    void register_variables();

  };

} //End namespace block

#endif // GRINS_INC_NAVIER_STOKES_BASE_H
