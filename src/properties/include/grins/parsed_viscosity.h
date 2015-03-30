//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_PARSED_VISCOSITY_H
#define GRINS_PARSED_VISCOSITY_H

//GRINS
#include "grins/viscosity_base.h"
#include "grins/assembly_context.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/function_base.h"

#include "libmesh/fem_system.h"

class GetPot;

namespace GRINS
{
  class ParsedViscosity : public ParameterUser,
                          public ViscosityBase
  {
  public:

    //! Deprecated constructor
    ParsedViscosity( const GetPot& input );
    virtual ~ParsedViscosity();

    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

    libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time=0 );

    void init(libMesh::FEMSystem* /*system*/){};

  private:

    ParsedViscosity();

    //! Helper function to ensure parsed function is non-zero
    void check_mu_nonzero( const std::string& function ) const;

    // User specified parsed function
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > _mu;

  };

  /* ------------------------- Inline Functions -------------------------*/  
  inline
    libMesh::Real ParsedViscosity::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // FIXME: We should be getting the variable index to get the qps from the context
    // not hardcode it to be 0
    const std::vector<libMesh::Point>& x = context.get_element_fe(0)->get_xyz();

    const libMesh::Point& x_qp = x[qp];

    libMesh::Number mu_value = (*_mu)(x_qp,context.time);

    return mu_value;
  }

  inline
  libMesh::Real ParsedViscosity::operator()( const libMesh::Point& p, const libMesh::Real time )
  {
    return (*_mu)(p,time);
  }

} // end namespace GRINS

#endif // GRINS_CONSTANT_VISCOSITY_H
