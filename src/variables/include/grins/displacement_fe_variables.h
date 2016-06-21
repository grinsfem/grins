//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_DISPLACEMENT_FE_VARIABLES_H
#define GRINS_DISPLACEMENT_FE_VARIABLES_H

// GRINS
#include "grins/multi_var_single_fe_type_variable.h"
#include "grins/variables_parsing.h"

namespace GRINS
{
  class DisplacementFEVariables : public MultiVarSingleFETypeVariable
  {
  public:

    //! Constructor
    /*! The arguments specify whether the spatial mesh is really 2D or 3D.
     *  This is needed for cases such as a 1D beam in 2D (is_2D = true)
     *  or 3D (is_3D = true) space or 2D shell manifolds in 3D (is_3D = true). */
    DisplacementFEVariables( const GetPot& input,
                             const std::string& physics_name,
                             bool is_2D, bool is_3D,
                             bool is_constraint_var = false );

    DisplacementFEVariables( const std::vector<std::string>& var_names,
                             const std::vector<VariableIndex>& var_indices )
      : MultiVarSingleFETypeVariable(var_names,var_indices),
        _have_v(false),
        _have_w(false),
        _u_idx(0),
        _v_idx(1),
        _w_idx(2),
        _is_2D(false),
        _is_3D(false)
    {
      if( var_names.size() > 1 )
        {
          _is_2D = true;
          _have_v = true;
        }

      if( var_names.size() > 2 )
        {
          _is_3D = true;
          _have_w = true;
        }
    }

    virtual ~DisplacementFEVariables(){};

    //! Initialize System variables
    virtual void init( libMesh::FEMSystem* system );

    bool have_v() const;
    bool have_w() const;

    VariableIndex u() const;
    VariableIndex v() const;
    VariableIndex w() const;

    const std::string& u_name() const;
    const std::string& v_name() const;
    const std::string& w_name() const;

  private:

    std::string subsection() const
    { return VariablesParsing::displacement_section(); }

    std::vector<std::string> old_var_names()
    {
      std::vector<std::string> var_names(3);
      var_names[0] = "u_displacment";
      var_names[1] = "v_displacment";
      var_names[2] = "w_displacment";
      return var_names;
    }

    std::vector<std::string> default_names()
    {
      std::vector<std::string> var_names(3);
      var_names[0] = "u";
      var_names[1] = "v";
      var_names[2] = "w";
      return var_names;
    }

    DisplacementFEVariables();

    bool _have_v;
    bool _have_w;

    unsigned int _u_idx, _v_idx, _w_idx;

    //! Tracks whether this is a 2D problem
    bool _is_2D;

    //! Tracks whether this is a 3D problem
    bool _is_3D;

  };

  inline
  VariableIndex DisplacementFEVariables::u() const
  {
    return this->_vars[_u_idx];
  }

  inline
  VariableIndex DisplacementFEVariables::v() const
  {
    libmesh_assert(_have_v);
    return this->_vars[_v_idx];
  }

  inline
  VariableIndex DisplacementFEVariables::w() const
  {
    libmesh_assert(_have_w);
    return this->_vars[_w_idx];
  }

  inline
  bool DisplacementFEVariables::have_v() const
  {
    return _have_v;
  }

  inline
  bool DisplacementFEVariables::have_w() const
  {
    return _have_w;
  }

  inline
  const std::string& DisplacementFEVariables::u_name() const
  {
    return this->_var_names[_u_idx];
  }

  inline
  const std::string& DisplacementFEVariables::v_name() const
  {
    return this->_var_names[_v_idx];
  }

  inline
  const std::string& DisplacementFEVariables::w_name() const
  {
    return this->_var_names[_w_idx];
  }

} // end namespace GRINS

#endif // GRINS_DISPLACEMENT_FE_VARIABLES_H
