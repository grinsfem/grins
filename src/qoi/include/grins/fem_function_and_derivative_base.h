//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_FEM_FUNCTION_AND_DERIVATIVE_BASE_H
#define GRINS_FEM_FUNCTION_AND_DERIVATIVE_BASE_H

// libMesh
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  /*!
    Extends libMesh::FEMFunctionBase by adding a function for computing derivatives
  */
  template<typename Output>
  class FEMFunctionAndDerivativeBase : public libMesh::FEMFunctionBase<Output>
  {
  public:

    //! Function Derivative Evaluation
    virtual void derivatives( libMesh::FEMContext & context,
                              const libMesh::Point & qp_xyz,
                              const libMesh::Real & JxW,
                              const unsigned int qoi_index,
                              const libMesh::Real time = 0.) = 0;

  };

}
#endif //GRINS_FEM_FUNCTION_AND_DERIVATIVE_BASE_H
