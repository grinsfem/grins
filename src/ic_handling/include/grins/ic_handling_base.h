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
#ifndef GRINS_IC_HANDLING_BASE_H
#define GRINS_IC_HANDLING_BASE_H

//GRINS
#include "grins/variable_name_defaults.h"
#include "grins/var_typedefs.h"
#include "grins/cached_values.h"

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"
#include "libmesh/auto_ptr.h"

// libMesh forward declarations
namespace libMesh
{
  template <typename Scalar>
  class CompositeFunction;

  class FEMSystem;
}

namespace GRINS
{
  //! Base class for reading and handling initial conditions for physics classes
  class ICHandlingBase
  {
  public:

    ICHandlingBase(const std::string& physics_name);

    virtual ~ICHandlingBase();

    void attach_initial_func( const libMesh::FunctionBase<libMesh::Number>& initial_val );

    virtual void read_ic_data( const GetPot& input, const std::string& id_str,
                               const std::string& ic_str,
                               const std::string& var_str,
                               const std::string& value_str );

    //! Override this method to initialize any system-dependent data.
    /*! Override this method to, for example, cache a System variable
      number. */
    virtual void init_ic_data( const libMesh::FEMSystem& system,
                               libMesh::CompositeFunction<libMesh::Number>& all_ics );

    // User will need to implement these functions for IC handling
    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_ic_types( const libMesh::subdomain_id_type ic_id,
                                const std::string& ic_id_string,
                                const int ic_type,
                                const std::string& ic_vars_string,
                                const std::string& ic_value_string,
                                const GetPot& input );

    libMesh::FunctionBase<libMesh::Number>* get_ic_func() const;

  protected:

    std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > _ic_func;

    std::string _physics_name;

    enum IC_BASE{ PARSED = -2,
                  CONSTANT };

    std::vector<std::string> _subfunction_variables;
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::FunctionBase<libMesh::Number>*
  ICHandlingBase::get_ic_func() const
  {
    return _ic_func.get();
  }
}
#endif // GRINS_BC_HANDLING_BASE
