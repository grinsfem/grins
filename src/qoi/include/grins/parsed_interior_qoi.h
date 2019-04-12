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
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/fem_function_base.h"
#include "libmesh/auto_ptr.h"

namespace GRINS
{
  //! Parsed Interior QoI
  /*!
    This class implements a QoI that is an arbitrary integral of a
    parsed function on the interior of the domain.
  */
  class ParsedInteriorQoI : public QoIBase
  {
  public:

    //! Constructor
    /*! Constructor takes GetPot object to read any input options associated
      with this QoI */
    ParsedInteriorQoI( const std::string& qoi_name );

    virtual ~ParsedInteriorQoI();

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    //! Initialize local variables
    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

    virtual void init_context( AssemblyContext& context );

    //! Compute the qoi value.
    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index );

    //! Compute the qoi derivative with respect to the solution.
    virtual void element_qoi_derivative( AssemblyContext& context,
                                         const unsigned int qoi_index );

  protected:

    std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >
    qoi_functional;

    //! Manual copy constructor due to the UniquePtr
    ParsedInteriorQoI(const ParsedInteriorQoI& original);

  private:
    //! User never call default constructor.
    ParsedInteriorQoI();

  };

  inline
  bool ParsedInteriorQoI::assemble_on_interior() const
  {
    return true;
  }

  inline
  bool ParsedInteriorQoI::assemble_on_sides() const
  {
    return false;
  }
}
#endif //GRINS_PARSED_INTERIOR_QOI_H
