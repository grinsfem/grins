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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef GRINS_STAB_HELPER_H
#define GRINS_STAB_HELPER_H

// libMesh
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/fe_base.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMContext;
}

namespace GRINS
{
  class StabilizationHelper
  {
  public:
    StabilizationHelper();
    ~StabilizationHelper();

    /*! \todo Should we inline this? */
    libMesh::RealGradient compute_g( libMesh::FEBase* fe,
				     libMesh::FEMContext& c,
				     unsigned int qp ) const;
    
    /*! \todo Should we inline this? */
    libMesh::RealTensor compute_G( libMesh::FEBase* fe,
				   libMesh::FEMContext& c,
				   unsigned int qp ) const;

  };

} // namespace GRINS
#endif // GRINS_STAB_HELPER_H
