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


#ifndef GRINS_INC_NAVIER_STOKES_H
#define GRINS_INC_NAVIER_STOKES_H

//GRINS
#include "grins/inc_navier_stokes_base.h"
#include "grins/pressure_pinning.h"

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
  class IncompressibleNavierStokes : public IncompressibleNavierStokesBase<Viscosity>
  {
  public:

    IncompressibleNavierStokes(const std::string& physics_name, const GetPot& input);

    ~IncompressibleNavierStokes(){};

    virtual void auxiliary_init( MultiphysicsSystem& system );

    //! Register postprocessing variables for IncompressibleNavierStokes
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    // Constraint part(s)
    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    //! Compute value of postprocessed quantities at libMesh::Point.
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

  protected:

    PressurePinning _p_pinning;

    //! Enable pressure pinning
    bool _pin_pressure;

    //! Index from registering this quantity
    unsigned int _mu_index;

  private:
    IncompressibleNavierStokes();

  };

} //End namespace block

#endif // GRINS_INC_NAVIER_STOKES_H
