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


#ifndef GRINS_PARSED_INTERIOR_QOI_H
#define GRINS_PARSED_INTERIOR_QOI_H

// GRINS
#include "grins/parsed_qoi_base.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  //! Parsed Interior QoI
  /*!
    This class implements a QoI that is an arbitrary integral of a
    parsed function on the interior of the domain.
  */
  class ParsedInteriorQoI : public ParsedQoIBase
  {
  public:

    using ParsedQoIBase::ParsedQoIBase;

    virtual ~ParsedInteriorQoI() = default;

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual QoIBase* clone() const override;

    virtual bool assemble_on_interior() const override
    { return true; }

    virtual bool assemble_on_sides() const override
    {  return false; }

    //! Initialize local variables
    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num ) override;

    virtual void init_context( AssemblyContext& context ) override;

    virtual void register_active_vars( std::set<unsigned int> & element_vars,
                                       std::set<unsigned int> & side_vars ) override;

    //! Compute the qoi value.
    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index ) override;

    //! Compute the qoi derivative with respect to the solution.
    virtual void element_qoi_derivative( AssemblyContext& context,
                                         const unsigned int qoi_index ) override;

  protected:

    //! We can use the default copy constructor, but still want it protected
    ParsedInteriorQoI(const ParsedInteriorQoI& original) = default;

  };
}
#endif //GRINS_PARSED_INTERIOR_QOI_H
