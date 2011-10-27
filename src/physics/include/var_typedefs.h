//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef VAR_TYPEDEFS_H
#define VAR_TYPEDEFS_H

#include <string>
#include <map>

namespace GRINS
{
  //! More descriptive name of the type used for (owned) variable indices
  typedef unsigned int VariableIndex;

  //! More descriptive name of the type used for (registered) variable indices
  typedef unsigned int RegtdVariableIndex;

  //! Map between variable name and system index
  typedef std::map<std::string,VariableIndex> VariableMap;
  typedef std::map<std::string,VariableIndex>::iterator VariableMapIt;
  typedef std::map<std::string,VariableIndex>::const_iterator VariableMapConstIt;
}
#endif //VAR_TYPEDEFS_H
