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

#ifndef GRINS_DIRICHLET_BC_FACTORY_ABSTRACT_H
#define GRINS_DIRICHLET_BC_FACTORY_ABSTRACT_H

// GRINS
#include "grins/bc_factory_abstract.h"

// libMesh
#include "libmesh/dirichlet_boundaries.h"

namespace GRINS
{
  class DirichletBCFactoryAbstract : public BCFactoryAbstract<libMesh::DirichletBoundary>
  {
  public:
    DirichletBCFactoryAbstract( const std::string& bc_type_name )
      : BCFactoryAbstract<libMesh::DirichletBoundary>(bc_type_name)
    {}

    ~DirichletBCFactoryAbstract(){};

  protected:

    //! Helper function
    /*! This will search the given section to make sure at least one of the var names as been
      set. An error will be thrown if no vars were found to be set. Those vars that have been
      set in the input file are returned in the vars_found variable, if vars_found is not NULL. */
    void check_for_vars( const GetPot& input, const std::string& section,
                         const std::vector<std::string>& var_names,
                         std::set<std::string>* vars_found = NULL );

  };

} // end namespace GRINS

#endif // GRINS_DIRICHLET_BC_FACTORY_ABSTRACT_H
