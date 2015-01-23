//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_SPALART_ALLMARAS_H
#define GRINS_SPALART_ALLMARAS_H

//GRINS
#include "grins/physics.h"
#include "primitive_flow_fe_variables.h"
#include "grins/turbulence_fe_variables.h"
#include "grins/turbulence_models_base.h"
#include "grins/spalart_allmaras_helper.h"

//Utils
#include "grins/distance_function.h"

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
  class SpalartAllmaras : public TurbulenceModelsBase<Viscosity>
  {
  public:

    SpalartAllmaras(const std::string& physics_name, const GetPot& input);

    ~SpalartAllmaras();
        
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );    

    // Element time derivative
    virtual void element_time_derivative(bool compute_jacobian, AssemblyContext& context, CachedValues& /*cache*/);
        
    // A distance function to get distances from boundaries to qps
    libMesh::AutoPtr<DistanceFunction> distance_function;

  protected:
    
    //! Physical dimension of problem
    /*! \todo Do we really need to cache this? */
    unsigned int _dim;

    // The flow variables
    PrimitiveFlowFEVariables _flow_vars;

    // These are defined for each physics
    TurbulenceFEVariables _turbulence_vars;   

    // Spalart Allmaras Helper object
    SpalartAllmarasHelper _spalart_allmaras_helper;

  private:
    SpalartAllmaras();

  };

} //End namespace block

#endif // GRINS_SPALART_ALLMARAS_H
