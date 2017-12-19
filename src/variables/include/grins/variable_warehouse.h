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


#ifndef GRINS_VARIABLE_WAREHOUSE_H
#define GRINS_VARIABLE_WAREHOUSE_H

// GRINS
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/libmesh_common.h"

namespace GRINS
{
  //! Not intended for users as a public API
  namespace GRINSPrivate
  {
    //! Track what FEVariablesBase objects have been created
    /*! Several modules need to interact with the Variables
      in use. So, this object creates a place to register
      a Variable class.

      \todo Currently, we only allow only unique variable types. This will
      change in the future once we allow per-subdomain variables.
      Hence, we have this in GRINSPrivate since the API may
      change in the future. */
    class VariableWarehouse
    {
    public:
      VariableWarehouse(){};

      ~VariableWarehouse(){};

      //! Check if variable is registered
      static bool is_registered( const std::string& var_name );

      //! First check if var_name is registered and then register
      /*! Use this API if you may be attempting to register the same
        variable more than once. */
      static void check_and_register_variable( const std::string& var_name,
                                               std::shared_ptr<FEVariablesBase>& variable );

      static void register_variable( const std::string& var_name,
                                     std::shared_ptr<FEVariablesBase>& variable );

      static std::shared_ptr<FEVariablesBase> get_variable_ptr( const std::string& var_name );

      static FEVariablesBase& get_variable( const std::string& var_name );

      template <typename DerivedType>
      static DerivedType& get_variable_subclass( const std::string& var_name );

      //! Clears the var_map()
      static void clear()
      { var_map().clear(); }

    protected:

      static std::map<std::string,std::shared_ptr<FEVariablesBase> >& var_map();

    };

    inline
    bool VariableWarehouse::is_registered( const std::string& var_name )
    {
      bool var_found = (var_map().find(var_name) != var_map().end() );
      return var_found;
    }

    inline
    void VariableWarehouse::check_and_register_variable( const std::string& var_name,
                                                         std::shared_ptr<FEVariablesBase>& variable )
    {
      if( !VariableWarehouse::is_registered(var_name) )
        VariableWarehouse::register_variable(var_name,variable);
    }

    inline
    void VariableWarehouse::register_variable( const std::string& var_name,
                                               std::shared_ptr<FEVariablesBase>& variable )
    {
      if( VariableWarehouse::is_registered(var_name) )
        libmesh_error_msg("ERROR: Duplicate FEVariable registration not allowed!");

      var_map()[var_name] = variable;
    }

    inline
    FEVariablesBase& VariableWarehouse::get_variable( const std::string& var_name )
    {
      std::shared_ptr<FEVariablesBase> var_ptr = VariableWarehouse::get_variable_ptr(var_name);
      return *var_ptr;
    }

    template <typename DerivedType>
    inline
    DerivedType& VariableWarehouse::get_variable_subclass( const std::string& var_name )
    {
      FEVariablesBase& var_base = VariableWarehouse::get_variable(var_name);

      DerivedType& derived_var = libMesh::cast_ref<DerivedType&>(var_base);

      return derived_var;
    }

  } // end namespace GRINSPrivate

} // end namespace GRINS

#endif // GRINS_VARIABLE_WAREHOUSE_H
