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

#ifndef GRINS_SPALART_ALLMARAS_STAB_BASE_H
#define GRINS_SPALART_ALLMARAS_STAB_BASE_H

//GRINS
#include "grins/spalart_allmaras.h"
#include "grins/spalart_allmaras_stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  template<class Viscosity>
  class SpalartAllmarasStabilizationBase : public SpalartAllmaras<Viscosity>
  {

  public:

    SpalartAllmarasStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~SpalartAllmarasStabilizationBase(){};

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    SpalartAllmarasStabilizationHelper _stab_helper;

  private:

    SpalartAllmarasStabilizationBase();

  }; // End SpalartAllmarasStabilizationBase class declarations

} // End namespace GRINS

#endif // GRINS_SPALART_ALLMARAS_STAB_BASE_H
