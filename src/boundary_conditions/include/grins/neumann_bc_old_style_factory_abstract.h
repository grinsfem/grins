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

#ifndef GRINS_NEUMANN_BC_OLD_STYLE_FACTORY_ABSTRACT_H
#define GRINS_NEUMANN_BC_OLD_STYLE_FACTORY_ABSTRACT_H

// GRINS
#include "grins/neumann_bc_factory_abstract.h"

namespace GRINS
{
  class NeumannBCOldStyleFactoryAbstract : public NeumannBCFactoryAbstract
  {
  public:
    NeumannBCOldStyleFactoryAbstract( const std::string& bc_type_name )
      : NeumannBCFactoryAbstract(bc_type_name)
    {}

    ~NeumannBCOldStyleFactoryAbstract(){};

    //! Input variable for parsing old style
    /*! Deprecated, only used for backward compatibility. */
    static void set_value_var_old_style( const std::string& value_var )
    { _value_var_old_style = value_var; }

    //! Input variable index for parsing old style
    /*! Deprecated, only used for backward compatibility. */
    static void set_value_index_old_style( unsigned int idx )
    { _value_idx_old_style = idx; }

  protected:

    //! Helper function to reduce code duplication
    virtual void check_state() const;

    //! Helper function to reduce code duplication
    virtual void reset_state();

    static std::string _value_var_old_style;

    static unsigned int _value_idx_old_style;

  };

} // end namespace GRINS

#endif // GRINS_NEUMANN_BC_OLD_STYLE_FACTORY_ABSTRACT_H
