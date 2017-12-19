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

// This class
#include "grins/variable_warehouse.h"

namespace GRINS
{
  namespace GRINSPrivate
  {
    std::map<std::string,std::shared_ptr<FEVariablesBase> >& VariableWarehouse::var_map()
    {
      static std::map<std::string,std::shared_ptr<FEVariablesBase> > _var_map;
      return _var_map;
    }

    std::shared_ptr<FEVariablesBase> VariableWarehouse::get_variable_ptr( const std::string& var_name )
    {
      if( !VariableWarehouse::is_registered(var_name) )
        {
          const std::map<std::string,std::shared_ptr<FEVariablesBase> >& map = var_map();

          std::stringstream error_msg;
          error_msg << "ERROR: Could not find Variable "+var_name+" in the VariableWarehouse!"
                    << std::endl
                    << "       Variables currently the VariableWarehouse are: "
                    << (map.begin())->first << std::endl;

          std::map<std::string,std::shared_ptr<FEVariablesBase> >::const_iterator it = map.begin();
          it++;
          for( ; it != map.end(); ++it )
            error_msg << std::string(54,' ') << it->first << std::endl;

          libmesh_error_msg(error_msg.str());
        }

      std::shared_ptr<FEVariablesBase> var_ptr = var_map()[var_name];

      if( !var_ptr )
        libmesh_error_msg("ERROR: Variable "+var_name+" is an invalid pointer!");

      return var_ptr;
    }

  } // end namespace GRINSPrivate
} // end namespace GRINS
