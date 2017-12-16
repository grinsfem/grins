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


#ifndef GRINS_SPALART_ALLMARAS_H
#define GRINS_SPALART_ALLMARAS_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/turbulence_models_base.h"
#include "grins/spalart_allmaras_helper.h"
#include "grins/spalart_allmaras_parameters.h"

//Utils
#include "grins/distance_function.h"

//libMesh
#include "libmesh/mesh.h"
#include "libmesh/boundary_info.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/boundary_mesh.h"

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

    ~SpalartAllmaras(){};

    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // Element time derivative
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    // A distance function to get distances from boundaries to qps
    std::unique_ptr<DistanceFunction> distance_function;

    // Boundary mesh objected that will be updated using the wall ids
    std::unique_ptr<libMesh::SerialMesh> boundary_mesh;

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    // The flow variables
    VelocityVariable& _flow_vars;
    PressureFEVariable& _press_var;
    // These are defined for each physics
    TurbulenceFEVariables& _turbulence_vars;

    // Spalart Allmaras Helper object
    SpalartAllmarasHelper _spalart_allmaras_helper;

    //! Object handling the plethora of parameters
    SpalartAllmarasParameters _sa_params;

    // Wall ids set, to be read in, tells us which bc_id's correspond to walls
    std::set<libMesh::boundary_id_type> _wall_ids;

    // No of walls
    unsigned int _no_of_walls;

    // Infinite distance case
    bool _infinite_distance;

  private:
    SpalartAllmaras();
  };

} //End namespace block

#endif // GRINS_SPALART_ALLMARAS_H
