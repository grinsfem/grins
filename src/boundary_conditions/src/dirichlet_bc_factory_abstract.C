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

// This class
#include "grins/dirichlet_bc_factory_abstract.h"

namespace GRINS
{
  // Full specialization for the Factory<DirichletBoundary>
  template<>
  std::map<std::string, FactoryAbstract<libMesh::DirichletBoundary>*>&
  FactoryAbstract<libMesh::DirichletBoundary>::factory_map()
  {
    static std::map<std::string, FactoryAbstract<libMesh::DirichletBoundary>*> _map;
    return _map;
  }

  // Definition of static members
  template<>
  const GetPot* FactoryWithGetPot<libMesh::DirichletBoundary>::_input = NULL;

  template<>
  MultiphysicsSystem* BCFactoryAbstract<libMesh::DirichletBoundary>::_system = NULL;

  template<>
  const std::set<BoundaryID>* BCFactoryAbstract<libMesh::DirichletBoundary>::_bc_ids = NULL;

  template<>
  const FEVariablesBase* BCFactoryAbstract<libMesh::DirichletBoundary>::_fe_var = NULL;

  template<>
  std::string BCFactoryAbstract<libMesh::DirichletBoundary>::_section = std::string("DIE!");

  void DirichletBCFactoryAbstract::check_for_vars( const GetPot& input, const std::string& section,
                                                   const std::vector<std::string>& var_names,
                                                   std::set<std::string>* vars_found )
  {
    if( vars_found )
      vars_found->clear();

    bool found_var = false;

    for( std::vector<std::string>::const_iterator name = var_names.begin();
         name < var_names.end(); ++name )
      {
        if( input.have_variable( section+"/"+(*name) ) )
          {
            found_var = true;

            if( vars_found )
              vars_found->insert( (*name) );
          }
      }

    // If no variables were found to be set, error out printing
    // the section and the variable names we were looking for.
    if( !found_var )
      {
        std::string err_msg = "ERROR: Could find any active variable assigned a Dirichlet boundary value\n";
        err_msg += "       in section "+section+". Active variables are:\n";

        for( std::vector<std::string>::const_iterator name = var_names.begin();
             name < var_names.end(); ++name )
          err_msg += "       "+(*name)+"\n";

        libmesh_error_msg(err_msg);
      }

    if( vars_found )
      libmesh_assert(!vars_found->empty());
  }

} // end namespace GRINS
