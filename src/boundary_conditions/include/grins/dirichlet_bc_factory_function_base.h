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

#ifndef GRINS_DIRICHLET_BC_FACTORY_FUNCTION_BASE_H
#define GRINS_DIRICHLET_BC_FACTORY_FUNCTION_BASE_H

// GRINS
#include "grins/dirichlet_bc_factory_abstract.h"
#include "grins/multiphysics_sys.h"

namespace GRINS
{
  template<typename FunctionType>
  class DirichletBCFactoryFunctionBase : public DirichletBCFactoryAbstract
  {
  public:

    DirichletBCFactoryFunctionBase( const std::string& bc_type_name )
      : DirichletBCFactoryAbstract(bc_type_name)
    {}

    ~DirichletBCFactoryFunctionBase(){};

  protected:

    //! Builds the FunctionBase object for boundary condition
    /*! Subclasses should override this function to build the
      FunctionBase object that corresponds to the variables passed
      in var_names. The variable names passed in will correspond to
      only a single VariableBase object, e.g.  Velocity. The section
      arguments corresponds to the section to parse for the
      variables in the input file, e.g. input(section+"/"+var_names[0]).

      The variable names passed in correspond to all the active
      variable names for a Variable group. The variable names
      present in var_names after this function is called will
      correspond to those variables. Note that var_names is *NOT*
      const. This is because one of the behaviors of subclasses may
      be to *remove* variables from the list. For example, for some
      symmetry conditions, we only want to enforce zero on certain
      components of the solution while leaving others untouched. */
    virtual std::unique_ptr<FunctionType>
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section ) =0;

    //! Dispatch, based on FunctionType, to the correct DirchletBoundary construction
    std::unique_ptr<libMesh::DirichletBoundary>
    make_dirichlet_boundary( const std::set<BoundaryID>& bc_ids,
                             const libMesh::System& system,
                             std::unique_ptr<FunctionType>& func,
                             const std::vector<VariableIndex>& var_indices );

    //! Helper function that can be overridded in subclasses
    /*! Mainly for backward compatibility with OldStyle parsing. */
    virtual const std::vector<std::string>& get_var_names() const
    { return this->_fe_var->active_var_names(); }

  private:

    virtual std::unique_ptr<libMesh::DirichletBoundary> create();

  };

} // end namespace GRINS

#endif // GRINS_DIRICHLET_BC_FACTORY_FUNCTION_BASE_H
