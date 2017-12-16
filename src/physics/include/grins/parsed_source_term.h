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

#ifndef GRINS_PARSED_SOURCE_TERM_H
#define GRINS_PARSED_SOURCE_TERM_H

// GRINS
#include "grins/source_term_base.h"

// libMesh
#include "libmesh/parsed_function.h"

namespace GRINS
{

  class ParsedSourceTerm : public SourceTermBase
  {
  public:

    ParsedSourceTerm( const std::string& physics_name, const GetPot& input );

    virtual ~ParsedSourceTerm(){}

    virtual void init_context( AssemblyContext& context );

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

  protected:

    libMesh::ParsedFunction<libMesh::Real> _value;

  private:

    ParsedSourceTerm();

  };

} // end namespace GRINS

#endif // GRINS_PARSED_SOURCE_TERM_H
