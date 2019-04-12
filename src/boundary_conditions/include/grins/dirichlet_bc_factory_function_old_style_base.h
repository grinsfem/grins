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

#ifndef GRINS_DIRICHLET_BC_FACTORY_FUNCTION_OLD_STYLE_BASE_H
#define GRINS_DIRICHLET_BC_FACTORY_FUNCTION_OLD_STYLE_BASE_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"

namespace GRINS
{
  template<typename FunctionType>
  class DirichletBCFactoryFunctionOldStyleBase : public DirichletBCFactoryFunctionBase<FunctionType>
  {
  public:

    DirichletBCFactoryFunctionOldStyleBase( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase<FunctionType>(bc_type_name)
    {}

    ~DirichletBCFactoryFunctionOldStyleBase(){};

    //! Input variable for parsing old style
    /*! Deprecated, only used for backward compatibility. */
    static void set_value_var_old_style( const std::string& value_var )
    { _value_var_old_style = value_var; }

    //! Input variable index for parsing old style
    /*! Deprecated, only used for backward compatibility. */
    static void set_value_index_old_style( unsigned int idx )
    { _value_idx_old_style = idx; }

    static void set_var_names_old_style( const std::vector<std::string>& var_names )
    { _var_names_old_style = &var_names; }

  protected:

    //! Helper function
    virtual void check_state() const;

    //! Helper function
    virtual void reset_state();

    virtual const std::vector<std::string>& get_var_names() const
    { return *(this->_var_names_old_style); }

    static std::string _value_var_old_style;

    static unsigned int _value_idx_old_style;

    static const std::vector<std::string>* _var_names_old_style;

  };

  template<typename FunctionType>
  inline
  void DirichletBCFactoryFunctionOldStyleBase<FunctionType>::check_state() const
  {
    DirichletBCFactoryAbstract::check_state();

    if( this->_value_var_old_style == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_value_var_old_style() before building boundary condition!");

    if( this->_value_idx_old_style == libMesh::invalid_uint )
      libmesh_error_msg("ERROR: must call set_value_index_old_style() before building boundary condition!");

    if( !this->_var_names_old_style )
      libmesh_error_msg("ERROR: must call set_var_names_old_style() before building boundary condition!");
  }

  template<typename FunctionType>
  inline
  void DirichletBCFactoryFunctionOldStyleBase<FunctionType>::reset_state()
  {
    DirichletBCFactoryAbstract::reset_state();
    this->_value_var_old_style = std::string("DIE!");
    this->_value_idx_old_style = libMesh::invalid_uint;
    this->_var_names_old_style = NULL;
  }

} // end namespace GRINS

#endif // GRINS_DIRICHLET_BC_FACTORY_FUNCTION_OLD_STYLE_BASE_H
