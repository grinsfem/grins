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


#ifndef GRINS_PARSED_QOI_BASE_H
#define GRINS_PARSED_QOI_BASE_H

// GRINS
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/fem_function_base.h"

// C++
#include <set>
#include <memory>

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  class ParsedQoIBase : public QoIBase
  {
  public:

    using QoIBase::QoIBase;

    virtual ~ParsedQoIBase() = default;

  protected:

    std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >
    qoi_functional;

    void init_qoi_functional( const GetPot & input,
                              const MultiphysicsSystem & system,
                              const std::string & input_string );

    //! Manual copy constructor due to the unique_ptr
    ParsedQoIBase(const ParsedQoIBase & original);
  };

} // end namespace GRINS

#endif // GRINS_PARSED_QOI_BASE_H
