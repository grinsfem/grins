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


#ifndef GRINS_TURBULENCE_MODELS_BASE_H
#define GRINS_TURBULENCE_MODELS_BASE_H

//GRINS
#include "grins/physics.h"
#include "grins/single_variable.h"

//Utils
#include "grins/distance_function.h"

namespace GRINS
{

  //! Physics class for Turbulence Models
  /*!
    This physics class implements the various turbulence models that can be used alongside
    NS to model turbulent flows.
    This is a templated class, the class Viscosity can be instantiated as a specific type
    (right now:ConstantViscosity or SpatiallyVaryingViscosity) to allow the user
    to specify a constant or spatially varying viscosity in the input file
  */
  template<class Viscosity>
  class TurbulenceModelsBase : public Physics
  {
  public:

    TurbulenceModelsBase(const std::string& physics_name, const GetPot& input);

    ~TurbulenceModelsBase(){};

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    //! Physical dimension of problem
    /*! \todo Do we really need to cache this? */
    unsigned int _dim;

    //! Material parameters, read from input
    /** \todo Create objects to allow for function specification */
    libMesh::Number _rho;

    //! Viscosity object
    Viscosity _mu;

  private:
    TurbulenceModelsBase();

  };

} //End namespace block

#endif // GRINS_TURBULENCE_MODELS_BASE_H
