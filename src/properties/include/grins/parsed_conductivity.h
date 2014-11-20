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


#ifndef GRINS_PARSED_CONDUCTIVITY_H
#define GRINS_PARSED_CONDUCTIVITY_H

//GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/function_base.h"

class GetPot;

namespace GRINS
{
  class ParsedConductivity
  {
  public:

    ParsedConductivity( const GetPot& input );
    ~ParsedConductivity();
    
    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;
    
  private:

     ParsedConductivity();
    
    // User specified parsed function
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > k;

  };

  /* ------------------------- Inline Functions -------------------------*/  
  inline
    libMesh::Real ParsedConductivity::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // FIXME: We should be getting the variable index to get the qps from the context
    // not hardcode it to be 0
    const std::vector<libMesh::Point>& x = context.get_element_fe(0)->get_xyz();

    const libMesh::Point& x_qp = x[qp];

    libMesh::Number _k_value = (*k)(x_qp,context.time);

    return _k_value;
  }    

} // end namespace GRINS

#endif // GRINS_PARSED_CONDUCTIVITY_H
