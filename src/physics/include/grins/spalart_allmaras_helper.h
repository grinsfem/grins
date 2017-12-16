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


#ifndef GRINS_SPALART_ALLMARAS_HELPER_H
#define GRINS_SPALART_ALLMARAS_HELPER_H

//GRINS
#include "grins/physics.h"

//Utils
#include "grins/distance_function.h"

namespace GRINS
{
  // Forward declarations
  class VelocityVariable;
  class PressureFEVariable;

  class SpalartAllmarasHelper
  {
  public:

    SpalartAllmarasHelper(const GetPot& input);

    virtual ~SpalartAllmarasHelper(){};

    void init_variables( libMesh::FEMSystem* system );

    // The vorticity function
    libMesh::Real vorticity(AssemblyContext& context, unsigned int qp) const;

  protected:

    // The flow variables
    const VelocityVariable& _flow_vars;
    const PressureFEVariable& _press_var;

  private:

    SpalartAllmarasHelper();

  };

} //End namespace block

#endif // GRINS_SPALART_ALLMARAS_HELPER_H
