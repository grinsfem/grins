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
#include "grins/neumann_bc_old_style_factory_abstract.h"

namespace GRINS
{
  std::string NeumannBCOldStyleFactoryAbstract::_value_var_old_style = std::string("DIE!");

  unsigned int NeumannBCOldStyleFactoryAbstract::_value_idx_old_style = libMesh::invalid_uint;

  void NeumannBCOldStyleFactoryAbstract::check_state() const
  {
    NeumannBCFactoryAbstract::check_state();

    if( this->_value_var_old_style == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_value_var_old_style() before building boundary condition!");

    if( this->_value_idx_old_style == libMesh::invalid_uint )
      libmesh_error_msg("ERROR: must call set_value_index_old_style() before building boundary condition!");
  }

  void NeumannBCOldStyleFactoryAbstract::reset_state()
  {
    NeumannBCFactoryAbstract::reset_state();
    this->_value_var_old_style = std::string("DIE!");
    this->_value_idx_old_style = libMesh::invalid_uint;
  }

} // end namespace GRINS
