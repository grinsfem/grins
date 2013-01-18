//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef VARIABLE_NAME_DEFAULTS_H
#define VARIABLE_NAME_DEFAULTS_H

#include <string>

namespace GRINS
{
  //! Default physics variable names
  /*!
    These are the default string names for all the available physics
    variables. These can be reset by the user, but we provide sane defaults
    here.
   */
  /** \todo Should we put the default physics variable names in a 
            class instead of just GRINS namespace? */
  //! x-velocity
  const std::string u_var_name_default = "u";
  
  //! y-velocity
  const std::string v_var_name_default = "v";
  
  //! z-velocity
  const std::string w_var_name_default = "w";

  //! r-velocity (axisymmetric case)
  const std::string u_r_var_name_default = "r_vel";

  //! z-velocity (axisymmetric case)
  const std::string u_z_var_name_default = "z_vel";

  //! pressure
  const std::string p_var_name_default = "p";

  //! temperature
  const std::string T_var_name_default = "T";
  
  //! Ex field
  const std::string Ex_var_name_default = "Ex";
  
  //! Ey field
  const std::string Ey_var_name_default = "Ey";

  //! Ez field
  const std::string Ez_var_name_default = "Ez";

  //! Bx field
  const std::string Bx_var_name_default = "Bx";
  
  //! By field
  const std::string By_var_name_default = "By";

  //! Bz field
  const std::string Bz_var_name_default = "Bz";
}

#endif //VARIABLE_NAME_DEFAULTS_H
