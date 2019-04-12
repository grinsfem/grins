//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_ABSORPTION_COEFF_BASE_H
#define GRINS_ABSORPTION_COEFF_BASE_H

// libMesh
#include "libmesh/fem_function_base.h"
#include "libmesh/auto_ptr.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/hitran.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/fem_function_and_derivative_base.h"

namespace GRINS
{
  class AbsorptionCoeffBase : public FEMFunctionAndDerivativeBase<libMesh::Real>
  {
  public:
    AbsorptionCoeffBase(libMesh::Real nu)
      : _nu(nu)
    {}

    libMesh::Real get_wavenumber() const
    {
      return _nu;
    }

    void set_wavenumber(libMesh::Real new_nu)
    {
      _nu = new_nu;
    }

    AbsorptionCoeffBase() = delete;

  protected:
    //! Desired wavenumber [cm^-1]
    libMesh::Real _nu;

  };

}
#endif //GRINS_ABSORPTION_COEFF_BASE_H
