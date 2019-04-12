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

#ifndef GRINS_NEUMANN_BC_FUNCTION_BASE_H
#define GRINS_NEUMANN_BC_FUNCTION_BASE_H

// GRINS
#include "grins/neumann_bc_abstract.h"
#include "grins/var_typedefs.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/auto_ptr.h" // UniquePtr
#include "libmesh/point.h"
#include "libmesh/function_base.h"
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  template<typename FunctionType, typename FEShape>
  class NeumannBCFunctionBase : public NeumannBCAbstract
  {
  public:

    //! Constructor for function with only one variable
    NeumannBCFunctionBase( VariableIndex var )
      : NeumannBCAbstract(),
        _vars(1,var)
    {}

    //! Constructor for function with several variables
    /*! This intended for FunctionTypes that are Composite so we can
      treat "group" variables, e.g. Displacement, as the same time
      using a composite function. */
    NeumannBCFunctionBase( const std::vector<VariableIndex>& vars )
      : NeumannBCAbstract(),
        _vars(vars)
    {}

    virtual ~NeumannBCFunctionBase(){};

    virtual bool eval_flux( bool compute_jacobian,
                            AssemblyContext& context,
                            libMesh::Real sign,
                            bool is_axisymmetric );

  protected:

    //! Helper function to dispatch to FEMFunctionBase API
    FEShape eval_func( AssemblyContext& context, const libMesh::Point& point,
                       libMesh::Real time, unsigned int component,
                       libMesh::FEMFunctionBase<FEShape>& func )
    {
      return func.component(context,component,point,time);
    }

    //! Helper function to dispatch to FunctionBase API
    FEShape eval_func( AssemblyContext& /*context*/, const libMesh::Point& point,
                       libMesh::Real time, unsigned int component,
                       libMesh::FunctionBase<FEShape>& func )
    {
      return func.component(component,point,time);
    }

    //! Variable indices for the variables whose Neumann contribution we're computing
    /*! We're assuming all _vars use the same FunctionType object. This allows us to
      use Composite type functions and treat "group" variables, like Displacement,
      at the same time. */
    std::vector<VariableIndex> _vars;

    //! Function object for the actual Neumann flux
    /*! Subclasses should initialize this function appropriately at construction
      time. */
    std::unique_ptr<FunctionType> _func;

  };

} // end namespace GRINS

#endif // GRINS_NEUMANN_BC_FUNCTION_BASE_H
